# 204 — Multi-sheet fan-in vs closed-loop contradiction

## Why this lane
You pointed out the key model mismatch:
- Closed-loop theorem assumes one sender/one receiver cycle.
- Copy-on-write can create multiple independent future senders.
- If those senders can route to the same origin branch, accumulation is fan-in,
  not a single cycle.

## New probe
- `ctc_multisheet_fanin_probe`

Outputs:
- `/tmp/bh_renders/ctc_multisheet_fanin_probe/ctc_multisheet_fanin_probe.txt`
- `/tmp/bh_renders/ctc_multisheet_fanin_probe/ctc_multisheet_fanin_probe.json`

Model:
- Per-sender send capability:
  - `s_{k+1} = max(0, eta*infra_gain*s_k + base_inflow - loss)`
- Sheet count:
  - `N_k = branching^k`
- Contribution to tracked origin:
  - split (no merge): `c_k = s_k`
  - fan-in (merge): `c_k = merge_fraction * N_k * s_k`

Conservation is enforced per generation:
- `c_k <= total_sent_k = N_k * s_k` for `merge_fraction <= 1`.

## Result
Default run (`generations=60`, `branching=2`, `merge_fraction=1`,
`eta=0.98`, `infra_gain=1.02`, `base_inflow=1e-9`, `loss=1e-10`):
- `per_sender_multiplier = eta*infra_gain = 0.9996` (near-neutral per sender)
- `fanin_seed_multiplier = branching*merge_fraction*eta*infra_gain = 1.9992`
- cumulative origin:
  - split: `1.694e-6 J`
  - fan-in: `1.233e11 J`

So your claim is correct in the model:
- many independent conservative sheets can produce huge accumulation at one
  receiver branch if many-to-one routing exists.

## Consistency with previous theorem
No contradiction with finding 202/203:
- Closed-cycle theorem blocks positive export for one closed no-drawdown cycle.
- Fan-in lane is not that model; it is multiple independent open contributors.

## New criterion
For fan-in seed growth, a useful indicator is:
- `R_fanin = branching * merge_fraction * eta * infra_gain`

Heuristic regimes:
- `R_fanin < 1`: fan-in seed contribution decays.
- `R_fanin = 1`: marginal.
- `R_fanin > 1`: fan-in growth strengthens (until capped by external limits).

Hard boundary remains model-dependent:
- Whether many-to-one routing to the same origin branch is physically allowed
  in the specific CTC realization.
