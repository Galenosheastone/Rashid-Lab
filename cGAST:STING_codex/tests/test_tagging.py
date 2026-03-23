from pathlib import Path

from mtgdeckdash.models import CardData, DeckCard
from mtgdeckdash.tags import TaggingEngine, apply_tags


def test_default_rules_apply_expected_tags() -> None:
    rules_path = Path(__file__).resolve().parents[1] / "mtgdeckdash" / "rules" / "default_rules.yaml"
    engine = TaggingEngine.from_yaml(rules_path)

    draw_card = CardData(
        name="Test Draw",
        oracle_text="Draw two cards.",
        type_line="Sorcery",
    )
    removal_card = CardData(
        name="Test Removal",
        oracle_text="Destroy target creature.",
        type_line="Instant",
    )

    assert "draw" in engine.tag_card(draw_card)
    assert "removal_single" in engine.tag_card(removal_card)


def test_apply_tags_on_deck_cards() -> None:
    rules_path = Path(__file__).resolve().parents[1] / "mtgdeckdash" / "rules" / "default_rules.yaml"
    engine = TaggingEngine.from_yaml(rules_path)

    deck_cards = [
        DeckCard(
            quantity=1,
            input_name="Counterspell",
            card=CardData(
                name="Counterspell",
                type_line="Instant",
                oracle_text="Counter target spell.",
            ),
        ),
        DeckCard(
            quantity=1,
            input_name="Island",
            card=CardData(name="Island", type_line="Basic Land — Island", oracle_text=""),
        ),
    ]

    tagged = apply_tags(deck_cards, engine)

    assert "counterspell" in tagged[0].tags
    assert tagged[1].tags == []
