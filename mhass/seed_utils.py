#!/usr/bin/env python3

import hashlib
PBSIM_SEED_MAX = 2_147_483_647

def derive_seed(master_seed, component, identifier=None):
    """
    Derive a deterministic 64-bit seed from a master seed and component name.

    Optionally include an identifier to derive object-specific seeds,
    such as one seed per PBSIM template.
    """
    parts = [str(master_seed), str(component)]

    if identifier is not None:
        parts.append(str(identifier))

    seed_material = "\x1f".join(parts).encode("utf-8")

    digest = hashlib.sha256(seed_material).digest()

    return int.from_bytes(digest[:8], byteorder="big", signed=False)


def derive_pbsim_seed(master_seed, template_identifier):
    """
    Derive a deterministic PBSIM-compatible seed for one template.

    The returned value is constrained to the positive signed 32-bit
    integer range: 1 through 2,147,483,647.
    """
    generic_seed = derive_seed(
        master_seed,
        "pbsim",
        template_identifier,
    )

    return (generic_seed % PBSIM_SEED_MAX) + 1
