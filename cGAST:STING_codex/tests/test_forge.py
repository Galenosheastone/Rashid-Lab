from __future__ import annotations

from pathlib import Path

import pytest

from mtgdeckdash.forge import generate_forge_dck
from mtgdeckdash.models import DeckCount, ResolvedDeck


def test_generate_forge_dck_has_expected_sections_and_counts(tmp_path: Path) -> None:
    resolved = ResolvedDeck(
        deck_id="abc123",
        name="Test Commander",
        commander_name="Test Commander",
        cards=[
            DeckCount(name="Test Commander", quantity=1),
            DeckCount(name="Wastes", quantity=99),
        ],
        analysis_dir=tmp_path,
        tagged_cards=[],
    )

    dck_path = generate_forge_dck(resolved)
    text = dck_path.read_text(encoding="utf-8")

    assert "[metadata]" in text
    assert "Deck Type=Commander" in text
    assert "[commander]" in text
    assert "1 Test Commander" in text
    assert "[main]" in text
    assert "99 Wastes" in text


def test_generate_forge_dck_raises_on_invalid_main_count(tmp_path: Path) -> None:
    resolved = ResolvedDeck(
        deck_id="bad123",
        name="Test Commander",
        commander_name="Test Commander",
        cards=[
            DeckCount(name="Test Commander", quantity=1),
            DeckCount(name="Wastes", quantity=98),
        ],
        analysis_dir=tmp_path,
        tagged_cards=[],
    )

    with pytest.raises(ValueError):
        generate_forge_dck(resolved)
