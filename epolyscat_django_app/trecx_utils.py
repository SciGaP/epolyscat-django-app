import logging
import re
from typing import NamedTuple

logger = logging.getLogger(__name__)


class Item(NamedTuple):
    inputValue: str
    default: str
    docu: str

    @property
    def empty(self):
        """Return True if input file had no value for this item."""
        return self.inputValue in ("EMPTY", "NOT_FOUND")


class Linp:
    """Read-only in memory representation of linp file"""

    LINE_REGEX = re.compile(
        r"^([^\s]+):([^[]*)\[(\d+)\]=(.*?)(?: \\inputValue{(.*)})? \\default{(.*)} \\docu{(.*)} \\allow"
    )

    def __init__(self, linp_file) -> None:
        """linp_file can be a file path str or file object"""
        self.categories = {}
        if isinstance(linp_file, str):
            with open(linp_file, "r") as lines:
                self.read_lines(lines)
        else:
            self.read_lines(linp_file)

    def read_lines(self, lines):
        for line in lines:
            if isinstance(line, bytes):
                line = line.decode()
            m = self.LINE_REGEX.match(line)
            if m:
                (
                    line_category,
                    line_name,
                    line_num,
                    deprecated_input_value,
                    input_value,
                    default_value,
                    docu,
                ) = m.groups()
                category = self.categories.setdefault(line_category, {})
                name = category.setdefault(line_name, {})
                name[line_num] = Item(
                    input_value or deprecated_input_value, default_value, docu
                )

    def item(self, category: str, name: str, line: int) -> Item:
        """Return the linp Item or None if not found."""
        names = self.categories.get(category, {})
        lines = names.get(name, {})
        return lines.get(str(line), None)


def is_empty(data_file):
    """Return False if data_file has any non-comment, non-empty lines."""

    line = data_file.readline()
    while line:
        if isinstance(line, bytes):
            line = line.decode()
        if not line.startswith("#") and line.strip() != "":
            return False
        line = data_file.readline()
    return True
