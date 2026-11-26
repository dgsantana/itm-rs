# NTIA ITM reference (test-only)

NTIA's C++ Irregular Terrain Model, copied unchanged from https://github.com/NTIA/itm at commit 183ad95bd813a8be11009df396e1c631356864b2
(release v1.4 plus 43 commits; `src/` here as `cpp/`, and `include/`). `LICENSE.md` is from commit
62168e686e8af4a234cf3f8067d80b543a5d52ef, the latest when it was copied, which corrects the statute it cites to Title 17;
the sources are identical at both commits. Works of NTIA employees are not subject to copyright protection in the United
States (Title 17 USC 105); NTIA's disclaimer is in `LICENSE.md`. This crate is a dev-dependency only: the
library's tests compare the Rust port with it. It is never built into the library.
