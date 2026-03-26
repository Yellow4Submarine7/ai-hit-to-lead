# Decision 0001: Public WWP2-C2 Workflow

## Status

Accepted.

## Decision

The public repository standardizes on the WWP2-C2 9M workflow:

`ZINC22 -> Uni-Dock shortlist -> GNINA/RF rescoring -> active learning`

## Rationale

- This is the production workflow we want external users to understand and reproduce.
- A single canonical public path keeps the stage layout, schemas, and tests consistent.
- Non-public historical experiments are intentionally left out of the runnable repository.

## Consequences

- The public stage map is `00` through `03`.
- Runnable examples, templates, and tests must align with the canonical WWP2-C2 path.
- Historical internal experiments are out of scope for this repository.
