from ensemblelab import Ensemble, generate
from ensemblelab.optimizers.mmff import MMFFOptimizer


def test_unoptimized_ensemble_show_contains_default_sections() -> None:
    ensemble = generate("CCO", n_confs=2)

    result = ensemble.show()

    assert "Ensemble" in result
    assert "Energy: uncomputed" in result
    assert "Conformers" in result
    assert "N/A" in result


def test_ensemble_show_can_combine_optional_sections() -> None:
    ensemble = generate("CCO", n_confs=2)

    result = ensemble.show(history=True, metadata=True, conformers=False)

    assert "Workflow History" in result
    assert "1. generation" in result
    assert "Metadata" in result
    assert "Conformers\n" not in result


def test_conformer_show_reports_stored_fields() -> None:
    conformer = generate("O", n_confs=1).conformers[0]

    result = conformer.show()

    assert "Conformer 0" in result
    assert "Energy          N/A" in result
    assert "Atoms           3" in result


def test_optimized_show_uses_relative_energy_and_optimizer_history() -> None:
    optimized = MMFFOptimizer().optimize(generate("CCO", n_confs=2))

    result = optimized.show(history=True)

    assert "Delta E (kcal/mol)" in result
    assert "Optimization: MMFF" in result
    assert "1. generation" in result
    assert "2. optimization" in result


def test_show_supports_legacy_functional_optimizer_history() -> None:
    generated = generate("O", n_confs=1)
    ensemble = Ensemble(
        smiles=generated.smiles,
        molecule=generated.molecule,
        conformers=generated.conformers,
        metadata={
            "optimization_history": [{"method": "MMFF", "max_steps": 500}],
        },
    )

    result = ensemble.show(history=True, conformers=False)

    assert "1. optimization" in result
    assert "Method: MMFF" in result
