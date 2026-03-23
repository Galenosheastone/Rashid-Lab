from mtgdeckdash.parser import parse_decklist_text


def test_parse_commander_header_and_card_lines() -> None:
    deck_text = """
// COMMANDER
1 Karn, Legacy Reforged (MAT) 49

// MAINBOARD
19 Wastes (FIN) 309
1 Sol Ring
"""

    parsed = parse_decklist_text(deck_text)

    assert len(parsed.entries) == 3
    assert parsed.commander_hint == "Karn, Legacy Reforged"

    commander_entry = parsed.entries[0]
    assert commander_entry.is_commander_hint is True
    assert commander_entry.quantity == 1
    assert commander_entry.set_code == "MAT"
    assert commander_entry.collector_number == "49"

    wastes = parsed.entries[1]
    assert wastes.name == "Wastes"
    assert wastes.quantity == 19
    assert wastes.set_code == "FIN"
    assert wastes.collector_number == "309"


def test_parse_ignores_non_card_lines_and_spaces() -> None:
    deck_text = """
// SIDEBOARD
This line should be ignored

  1   Arcane Signet   
2   Island
"""

    parsed = parse_decklist_text(deck_text)

    assert len(parsed.entries) == 2
    assert parsed.entries[0].name == "Arcane Signet"
    assert parsed.entries[0].quantity == 1
    assert parsed.entries[1].name == "Island"
    assert parsed.entries[1].quantity == 2
