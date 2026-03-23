from __future__ import annotations

import os
import json
from pathlib import Path

import pytest

from mtgdeckdash.models import CardData, DeckAnalysis, DeckCard
from mtgdeckdash.models import DeckInput
from mtgdeckdash.simulation import run_matchup_simulation


def _write_analysis_artifacts(analysis_dir: Path, commander_name: str) -> None:
    analysis_dir.mkdir(parents=True, exist_ok=True)
    analysis = DeckAnalysis(commander_name=commander_name, total_cards=100)
    (analysis_dir / "summary.json").write_text(
        json.dumps(analysis.model_dump(mode="json"), indent=2),
        encoding="utf-8",
    )

    tagged_cards = [
        DeckCard(
            quantity=1,
            input_name=commander_name,
            card=CardData(name=commander_name, type_line="Legendary Creature"),
        ),
        DeckCard(
            quantity=99,
            input_name="Wastes",
            card=CardData(name="Wastes", type_line="Basic Land — Wastes"),
        ),
    ]
    (analysis_dir / "tagged_cards.json").write_text(
        json.dumps([c.model_dump(mode="json") for c in tagged_cards], indent=2),
        encoding="utf-8",
    )


def test_simulation_end_to_end_optional_docker(tmp_path: Path) -> None:
    if not os.environ.get("DOCKER_AVAILABLE"):
        pytest.skip("Set DOCKER_AVAILABLE=1 to run docker-backed simulation test.")

    analysis_a = tmp_path / "analysis_a"
    analysis_b = tmp_path / "analysis_b"
    _write_analysis_artifacts(analysis_a, commander_name="Commander A")
    _write_analysis_artifacts(analysis_b, commander_name="Commander B")

    deck_a = DeckInput(source_type="analysis_output_dir", source_value=str(analysis_a), label="Deck A")
    deck_b = DeckInput(source_type="analysis_output_dir", source_value=str(analysis_b), label="Deck B")

    out_dir = tmp_path / "sim_out"
    result = run_matchup_simulation(
        deck_a_input=deck_a,
        deck_b_input=deck_b,
        games=2,
        seed=1,
        timeout_seconds=120,
        out_dir=out_dir,
        cache_root=tmp_path / "cache",
    )

    assert Path(result["report_html"]).exists()
    assert Path(result["summary_json"]).exists()


def test_simulation_local_backend_without_docker(tmp_path: Path) -> None:
    analysis_a = tmp_path / "analysis_a"
    analysis_b = tmp_path / "analysis_b"
    _write_analysis_artifacts(analysis_a, commander_name="Commander A")
    _write_analysis_artifacts(analysis_b, commander_name="Commander B")

    deck_a = DeckInput(source_type="analysis_output_dir", source_value=str(analysis_a), label="Deck A")
    deck_b = DeckInput(source_type="analysis_output_dir", source_value=str(analysis_b), label="Deck B")

    out_dir = tmp_path / "sim_local"
    result = run_matchup_simulation(
        deck_a_input=deck_a,
        deck_b_input=deck_b,
        games=25,
        seed=7,
        timeout_seconds=120,
        out_dir=out_dir,
        cache_root=tmp_path / "cache",
        backend="local",
    )

    summary_payload = json.loads(Path(result["summary_json"]).read_text(encoding="utf-8"))
    assert summary_payload["backend"] == "local_heuristic"
    assert summary_payload["result"]["games_recorded"] == 25
    assert Path(result["report_html"]).exists()
