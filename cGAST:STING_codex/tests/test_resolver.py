from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from mtgdeckdash.models import CardData, DeckAnalysis, DeckCard, DeckInput
from mtgdeckdash.resolver import resolve_deck_input


def _write_analysis_artifacts(analysis_dir: Path, commander_name: str = "Test Commander") -> None:
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
            tags=[],
        ),
        DeckCard(
            quantity=99,
            input_name="Wastes",
            card=CardData(name="Wastes", type_line="Basic Land — Wastes"),
            tags=["ramp"],
        ),
    ]

    (analysis_dir / "tagged_cards.json").write_text(
        json.dumps([c.model_dump(mode="json") for c in tagged_cards], indent=2),
        encoding="utf-8",
    )

    pd.DataFrame(
        [
            {"quantity": 1, "input_name": commander_name, "resolved_name": commander_name, "tags": ""},
            {"quantity": 99, "input_name": "Wastes", "resolved_name": "Wastes", "tags": "ramp"},
        ]
    ).to_csv(analysis_dir / "cards.csv", index=False)


def test_resolve_deck_input_decklist_text(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    cache_root = tmp_path / "cache"
    calls = {"count": 0}

    def fake_run_analysis_pipeline(
        deck_text: str,
        output_dir: Path,
        commander_override: str | None = None,
        rules_path: Path | None = None,
        deck_format: str = "commander",
        progress_callback=None,
    ):
        calls["count"] += 1
        _write_analysis_artifacts(output_dir, commander_name=commander_override or "Text Commander")
        return {"report_path": output_dir / "report.html"}

    monkeypatch.setattr("mtgdeckdash.resolver.run_analysis_pipeline", fake_run_analysis_pipeline)

    deck_input = DeckInput(
        source_type="decklist_text",
        source_value="// COMMANDER\n1 Text Commander\n99 Wastes",
        label="Deck A",
        commander_override="Text Commander",
    )

    resolved = resolve_deck_input(deck_input, cache_root=cache_root)

    assert resolved.commander_name == "Text Commander"
    assert resolved.analysis_dir.exists()
    assert resolved.deck_id
    assert calls["count"] == 1

    resolved_again = resolve_deck_input(deck_input, cache_root=cache_root)
    assert resolved_again.deck_id == resolved.deck_id
    assert calls["count"] == 1


def test_resolve_deck_input_decklist_file(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    cache_root = tmp_path / "cache"
    deck_file = tmp_path / "deck_b.txt"
    deck_file.write_text("// COMMANDER\n1 File Commander\n99 Wastes", encoding="utf-8")

    def fake_run_analysis_pipeline(
        deck_text: str,
        output_dir: Path,
        commander_override: str | None = None,
        rules_path: Path | None = None,
        deck_format: str = "commander",
        progress_callback=None,
    ):
        _write_analysis_artifacts(output_dir, commander_name="File Commander")
        return {"report_path": output_dir / "report.html"}

    monkeypatch.setattr("mtgdeckdash.resolver.run_analysis_pipeline", fake_run_analysis_pipeline)

    deck_input = DeckInput(
        source_type="decklist_file",
        source_value=str(deck_file),
        label="Deck B",
    )

    resolved = resolve_deck_input(deck_input, cache_root=cache_root)

    assert resolved.commander_name == "File Commander"
    assert any(card.name == "Wastes" and card.quantity == 99 for card in resolved.cards)


def test_resolve_deck_input_analysis_output_dir(tmp_path: Path) -> None:
    cache_root = tmp_path / "cache"
    analysis_dir = tmp_path / "analysis_a"
    _write_analysis_artifacts(analysis_dir, commander_name="Output Commander")

    deck_input = DeckInput(
        source_type="analysis_output_dir",
        source_value=str(analysis_dir),
        label="Deck A",
    )

    resolved = resolve_deck_input(deck_input, cache_root=cache_root)

    assert resolved.analysis_dir == analysis_dir.resolve()
    assert resolved.commander_name == "Output Commander"
    assert len(resolved.tagged_cards) == 2
