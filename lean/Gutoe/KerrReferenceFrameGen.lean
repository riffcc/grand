/- 
 * GUTOE — Lean Kerr Reference Frame Generator
 * Emits a low-res CSV intensity frame from Lean reference equations.
 * If output path ends with `.ppm`, emits a native false-color PPM image.
-/

import Gutoe.KerrReferenceFrame

open Gutoe.KerrReferenceFrame

def parseNatOr (s : String) (fallback : Nat) : Nat :=
  match s.toNat? with
  | some n => n
  | none => fallback

def pow10 (n : Nat) : Float :=
  (List.range n).foldl (fun acc _ => acc * 10.0) 1.0

def parseFloatSimple? (s : String) : Option Float :=
  match s.splitOn "." with
  | [i] =>
      i.toNat?.map Float.ofNat
  | [i, f] =>
      match i.toNat?, f.toNat? with
      | some iv, some fv =>
          some (Float.ofNat iv + Float.ofNat fv / pow10 f.length)
      | _, _ => none
  | _ => none

def parseFloatOr (s : String) (fallback : Float) : Float :=
  match parseFloatSimple? s with
  | some x => x
  | none => fallback

def rowCsv
    (tileW tileH y : Nat)
    (fov rObs r_s aStar thetaObs : Float)
    (fullW fullH x0 y0 : Nat) : String :=
  let vals := (List.range tileW).map (fun x =>
    s!"{referencePixelFloat fullW fullH (x0 + x) (y0 + y) fov rObs r_s aStar thetaObs}")
  String.intercalate "," vals

def rowFloats
    (tileW tileH y : Nat)
    (fov rObs r_s aStar thetaObs : Float)
    (fullW fullH x0 y0 : Nat) : List Float :=
  (List.range tileW).map (fun x =>
    referencePixelFloat fullW fullH (x0 + x) (y0 + y) fov rObs r_s aStar thetaObs)

def frameGridParallelIO
    (tileW tileH : Nat)
    (fov rObs r_s aStar thetaObs : Float)
    (fullW fullH x0 y0 : Nat) : IO (Array (List Float)) := do
  let tasks ← (List.range tileH).mapM (fun y =>
    IO.asTask (pure (rowFloats tileW tileH y fov rObs r_s aStar thetaObs fullW fullH x0 y0)))
  let rows ← tasks.mapM (fun t => IO.ofExcept t.get)
  pure rows.toArray

def frameCsvFromGrid (grid : Array (List Float)) : String :=
  let rows := grid.toList.map (fun row =>
    String.intercalate "," (row.map (fun v => s!"{v}")))
  String.intercalate "\n" rows

def clamp01 (x : Float) : Float :=
  if x < 0.0 then 0.0 else if x > 1.0 then 1.0 else x

def floatByte (x : Float) : Nat :=
  UInt64.toNat (Float.toUInt64 (Float.round (clamp01 x * 255.0)))

def toneMap (x exposure gamma blackLevel : Float) : Float :=
  let x0 := clamp01 x
  let lifted := clamp01 ((x0 - blackLevel) / max (1.0 - blackLevel) 1e-6)
  let exposed := clamp01 (exposure * lifted)
  clamp01 (exposed ^ gamma)

def opticalColor (tRaw : Float) : Nat × Nat × Nat :=
  -- Warm "blackbody-like" false color:
  -- dark -> deep red -> amber -> near-white.
  let t := Float.sqrt (clamp01 tRaw)
  let r := clamp01 (1.35 * t)
  let g := clamp01 (0.92 * t * t + 0.08 * t)
  let b := clamp01 (0.22 * t * t * t)
  (floatByte r, floatByte g, floatByte b)

def rowPpm
    (row : List Float)
    (exposure gamma blackLevel : Float) : String :=
  let vals := row.map (fun v =>
    let tv := toneMap (v / 255.0) exposure gamma blackLevel
    let (r, g, b) := opticalColor tv
    s!"{r} {g} {b}")
  String.intercalate " " vals

def framePpm
    (tileW tileH : Nat)
    (grid : Array (List Float))
    (exposure gamma blackLevel : Float) : String :=
  let rowsList := grid.toList
  let rows := rowsList.map (fun row => rowPpm row exposure gamma blackLevel)
  let header := s!"P3\n{tileW} {tileH}\n255\n"
  header ++ String.intercalate "\n" rows ++ "\n"

def main (args : List String) : IO UInt32 := do
  let outPath := args.getD 0 "/tmp/kerr_reference_128.csv"
  let w := parseNatOr (args.getD 1 "128") 128
  let h := parseNatOr (args.getD 2 "128") 128

  -- M87-like defaults for reference snapshots.
  -- Optional CLI args:
  --   3: fov
  --   4: rObs
  --   5: r_s
  --   6: aStar
  --   7: thetaObs (radians)
  let fov : Float := parseFloatOr (args.getD 3 "14.0") 14.0
  let rObs : Float := parseFloatOr (args.getD 4 "40.0") 40.0
  let r_s : Float := parseFloatOr (args.getD 5 "1.0") 1.0
  let aStar : Float := parseFloatOr (args.getD 6 "0.9") 0.9
  let thetaObs : Float := parseFloatOr (args.getD 7 "0.296705972839036") (3.141592653589793 * 17.0 / 180.0)
  -- Optional tile args for process-level tiling/stitch:
  --   8: full image width
  --   9: full image height
  --  10: tile x offset in full image
  --  11: tile y offset in full image
  --  12: exposure
  --  13: gamma
  --  14: black level
  let fullW : Nat := parseNatOr (args.getD 8 s!"{w}") w
  let fullH : Nat := parseNatOr (args.getD 9 s!"{h}") h
  let x0 : Nat := parseNatOr (args.getD 10 "0") 0
  let y0 : Nat := parseNatOr (args.getD 11 "0") 0
  let exposure : Float := parseFloatOr (args.getD 12 "1.15") 1.15
  let gamma : Float := parseFloatOr (args.getD 13 "1.55") 1.55
  let blackLevel : Float := parseFloatOr (args.getD 14 "0.10") 0.10

  -- Compute the full grid in parallel using Lean tasks.
  let grid ← frameGridParallelIO w h fov rObs r_s aStar thetaObs fullW fullH x0 y0

  if outPath.endsWith ".ppm" then
    let ppm := framePpm w h grid exposure gamma blackLevel
    IO.FS.writeFile outPath ppm
    IO.println s!"wrote Lean Kerr reference frame (PPM): {outPath} ({w}x{h})"
  else
    let csv := frameCsvFromGrid grid
    IO.FS.writeFile outPath csv
    IO.println s!"wrote Lean Kerr reference frame (CSV): {outPath} ({w}x{h})"
  pure 0
