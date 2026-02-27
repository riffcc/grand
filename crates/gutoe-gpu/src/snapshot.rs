use std::fs;
use std::io;
use std::path::Path;

const MAGIC: [u8; 8] = *b"GUTOESNP";
const VERSION_V1: u16 = 1;
const HEADER_LEN: usize = 8 + 2 + 2 + 8 + 8 + 8 + 8;

#[derive(Debug, Clone, PartialEq)]
pub struct UniverseSnapshot {
    pub tick: u64,
    pub seed: u64,
    pub sim_time: f64,
    pub payload: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SnapshotError {
    TooShort,
    BadMagic,
    UnsupportedVersion(u16),
    Truncated,
    HashMismatch { expected: u64, got: u64 },
}

impl std::fmt::Display for SnapshotError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TooShort => write!(f, "snapshot too short"),
            Self::BadMagic => write!(f, "invalid snapshot magic"),
            Self::UnsupportedVersion(v) => write!(f, "unsupported snapshot version: {v}"),
            Self::Truncated => write!(f, "snapshot payload truncated"),
            Self::HashMismatch { expected, got } => {
                write!(
                    f,
                    "snapshot hash mismatch: expected {expected:#x}, got {got:#x}"
                )
            }
        }
    }
}

impl std::error::Error for SnapshotError {}

#[inline]
fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325_u64;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

pub fn encode_snapshot(snapshot: &UniverseSnapshot) -> Vec<u8> {
    let payload_len = snapshot.payload.len() as u64;
    let payload_hash = fnv1a64(&snapshot.payload);
    let mut out = Vec::with_capacity(HEADER_LEN + snapshot.payload.len());
    out.extend_from_slice(&MAGIC);
    out.extend_from_slice(&VERSION_V1.to_le_bytes());
    out.extend_from_slice(&0_u16.to_le_bytes()); // reserved/flags
    out.extend_from_slice(&snapshot.tick.to_le_bytes());
    out.extend_from_slice(&snapshot.seed.to_le_bytes());
    out.extend_from_slice(&snapshot.sim_time.to_le_bytes());
    out.extend_from_slice(&payload_len.to_le_bytes());
    out.extend_from_slice(&payload_hash.to_le_bytes());
    out.extend_from_slice(&snapshot.payload);
    out
}

pub fn decode_snapshot(bytes: &[u8]) -> Result<UniverseSnapshot, SnapshotError> {
    if bytes.len() < HEADER_LEN {
        return Err(SnapshotError::TooShort);
    }
    if bytes[0..8] != MAGIC {
        return Err(SnapshotError::BadMagic);
    }
    let version = u16::from_le_bytes([bytes[8], bytes[9]]);
    if version != VERSION_V1 {
        return Err(SnapshotError::UnsupportedVersion(version));
    }

    let mut o = 12; // skip magic/version/flags
    let read_u64 = |buf: &[u8], off: &mut usize| -> u64 {
        let v = u64::from_le_bytes(buf[*off..*off + 8].try_into().expect("slice len"));
        *off += 8;
        v
    };
    let tick = read_u64(bytes, &mut o);
    let seed = read_u64(bytes, &mut o);
    let sim_time = f64::from_le_bytes(bytes[o..o + 8].try_into().expect("slice len"));
    o += 8;
    let payload_len = read_u64(bytes, &mut o) as usize;
    let expected_hash = read_u64(bytes, &mut o);
    let end = o.saturating_add(payload_len);
    if end > bytes.len() {
        return Err(SnapshotError::Truncated);
    }
    let payload = bytes[o..end].to_vec();
    let got_hash = fnv1a64(&payload);
    if got_hash != expected_hash {
        return Err(SnapshotError::HashMismatch {
            expected: expected_hash,
            got: got_hash,
        });
    }
    Ok(UniverseSnapshot {
        tick,
        seed,
        sim_time,
        payload,
    })
}

pub fn write_snapshot_file(path: &Path, snapshot: &UniverseSnapshot) -> io::Result<()> {
    fs::write(path, encode_snapshot(snapshot))
}

pub fn read_snapshot_file(path: &Path) -> Result<UniverseSnapshot, io::Error> {
    let bytes = fs::read(path)?;
    decode_snapshot(&bytes).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snapshot_roundtrip() {
        let s = UniverseSnapshot {
            tick: 42,
            seed: 1337,
            sim_time: 12.5,
            payload: vec![1, 2, 3, 4, 5, 9, 8],
        };
        let b = encode_snapshot(&s);
        let got = decode_snapshot(&b).expect("decode");
        assert_eq!(got, s);
    }

    #[test]
    fn snapshot_rejects_bad_magic() {
        let s = UniverseSnapshot {
            tick: 1,
            seed: 2,
            sim_time: 3.0,
            payload: vec![9, 9, 9],
        };
        let mut b = encode_snapshot(&s);
        b[0] ^= 0xFF;
        let err = decode_snapshot(&b).expect_err("should fail");
        assert!(matches!(err, SnapshotError::BadMagic));
    }

    #[test]
    fn snapshot_rejects_payload_corruption() {
        let s = UniverseSnapshot {
            tick: 99,
            seed: 7,
            sim_time: 1.25,
            payload: vec![11, 22, 33, 44],
        };
        let mut b = encode_snapshot(&s);
        let n = b.len();
        b[n - 1] ^= 0xAA;
        let err = decode_snapshot(&b).expect_err("should fail");
        assert!(matches!(err, SnapshotError::HashMismatch { .. }));
    }
}
