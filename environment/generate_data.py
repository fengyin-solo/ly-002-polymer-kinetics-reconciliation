#!/usr/bin/env python3
"""Build the fixed laboratory fixture used by the polymerization recon task."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
from pathlib import Path


GAS_CONSTANT = 8.31446261815324
TRUE_ACTIVATION_ENERGY = 92500.0
TRUE_PREFACTOR = 1.18e9
OBSERVATION_TIMES = (720.0, 1440.0, 2160.0, 2880.0)


def thermal_dose(steps, elapsed, activation_energy):
    """Integrate exp(-Ea/RT) over a piecewise-linear temperature program."""
    remaining = float(elapsed)
    total = 0.0
    for step in steps:
        if remaining <= 0.0:
            break
        duration = min(remaining, float(step["duration_s"]))
        fraction = duration / float(step["duration_s"])
        start = float(step["temp_start_k"])
        end = start + (float(step["temp_end_k"]) - start) * fraction
        # Ten-point Gauss-Legendre on each ramp is stable for these smooth histories.
        nodes = (
            -0.9739065285171717,
            -0.8650633666889845,
            -0.6794095682990244,
            -0.4333953941292472,
            -0.1488743389816312,
            0.1488743389816312,
            0.4333953941292472,
            0.6794095682990244,
            0.8650633666889845,
            0.9739065285171717,
        )
        weights = (
            0.0666713443086881,
            0.1494513491505806,
            0.2190863625159820,
            0.2692667193099963,
            0.2955242247147529,
            0.2955242247147529,
            0.2692667193099963,
            0.2190863625159820,
            0.1494513491505806,
            0.0666713443086881,
        )
        for node, weight in zip(nodes, weights):
            temperature = start + (end - start) * 0.5 * (node + 1.0)
            total += 0.5 * duration * weight * math.exp(
                -activation_energy / (GAS_CONSTANT * temperature)
            )
        remaining -= duration
    return total


def reaction_steps(index, rng):
    """Create plausible ramp, soak, and trim stages for one vessel run."""
    base = 338.0 + (index % 7) * 2.4 + rng.uniform(-2.0, 2.0)
    peak = 382.0 + (index * 3 % 17) * 1.1 + rng.uniform(-1.4, 1.4)
    hold = 366.0 + (index * 5 % 13) * 1.3 + rng.uniform(-1.2, 1.2)
    first = 780.0 + (index % 5) * 70.0 + rng.choice((-35.0, -15.0, 20.0, 45.0))
    second = 1020.0 + (index % 4) * 65.0 + rng.choice((-25.0, 10.0, 35.0))
    third = 1100.0 + (index % 3) * 90.0 + rng.choice((-30.0, 25.0))
    return [
        {
            "duration_s": round(first, 3),
            "temp_start_k": round(base, 3),
            "temp_end_k": round(peak, 3),
        },
        {
            "duration_s": round(second, 3),
            "temp_start_k": round(peak, 3),
            "temp_end_k": round(hold, 3),
        },
        {
            "duration_s": round(third, 3),
            "temp_start_k": round(hold, 3),
            "temp_end_k": round(hold + rng.choice((-2.5, 0.8, 3.1)), 3),
        },
    ]


def history_identifier(index):
    token = hashlib.sha256(f"ly-002-history-{index}".encode()).hexdigest()[:7].upper()
    return f"RX-{token[:3]}-{token[3:]}"


def bundle_identifier(index):
    token = hashlib.sha256(f"ly-002-sample-{index}".encode()).hexdigest()[:7].upper()
    return f"SPL-{token[:3]}-{token[3:]}"


def make_histories(rng, history_count=26):
    histories = []
    for index in range(history_count):
        histories.append(
            {
                "history_id": history_identifier(index),
                "steps": reaction_steps(index, rng),
                "catalyst_scale": round(0.72 + (index * 7 % 19) / 25.0, 6),
            }
        )
    return histories


def observed_conversion(history, elapsed):
    dose = thermal_dose(history["steps"], elapsed, TRUE_ACTIVATION_ENERGY)
    exponent = history["catalyst_scale"] * TRUE_PREFACTOR * dose
    return 1.0 - math.exp(-exponent)


def make_bundles(histories, rng, orphan_count=2):
    bundles = []
    true_order = list(range(len(histories)))
    rng.shuffle(true_order)
    for bundle_index, history_index in enumerate(true_order):
        history = histories[history_index]
        observations = []
        for time_index, elapsed in enumerate(OBSERVATION_TIMES):
            conversion = observed_conversion(history, elapsed)
            # One aliquot in four has a handling incident, not a model residual.
            if time_index == (bundle_index * 5 + 3) % 4:
                conversion = min(0.987, conversion + rng.choice((0.105, -0.092, 0.137)))
            else:
                conversion += rng.uniform(-0.00018, 0.00018)
            observations.append(
                {
                    "elapsed_s": elapsed,
                    "conversion": round(max(0.002, min(0.987, conversion)), 8),
                }
            )
        bundles.append(
            {
                "bundle_id": bundle_identifier(bundle_index),
                "observations": observations,
            }
        )

    # Two sample trays lost their reactor association entirely. Their readings
    # are plausible but do not share one kinetic signature.
    for extra_index in range(orphan_count):
        left = histories[(extra_index * 9 + 4) % len(histories)]
        right = histories[(extra_index * 9 + 15) % len(histories)]
        observations = []
        for time_index, elapsed in enumerate(OBSERVATION_TIMES):
            if time_index % 2 == 0:
                conversion = observed_conversion(left, elapsed)
            else:
                conversion = observed_conversion(right, elapsed)
            observations.append(
                {
                    "elapsed_s": elapsed,
                    "conversion": round(min(0.975, max(0.01, conversion)), 8),
                }
            )
        bundles.append(
            {
                "bundle_id": bundle_identifier(len(histories) + extra_index),
                "observations": observations,
            }
        )
    rng.shuffle(bundles)
    return bundles


def build_dataset(history_count=26):
    rng = random.Random(20261022)
    histories = make_histories(rng, history_count)
    bundles = make_bundles(histories, rng)
    return {
        "schema_version": 1,
        "reaction_histories": histories,
        "sample_bundles": bundles,
        "measurement_policy": {
            "activation_energy_bounds_j_mol": [40000.0, 150000.0],
            "ln_prefactor_bounds": [16.0, 28.0],
            "conversion_floor": 1e-9,
            "consensus_min_points": 3,
            "conversion_tolerance": 0.0015,
            "max_assignments": history_count,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--history-count", type=int, default=26)
    parser.add_argument("--filename", default="polymer_batches.json")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    target = args.out / args.filename
    target.write_text(
        json.dumps(
            build_dataset(args.history_count),
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(target)


if __name__ == "__main__":
    main()
