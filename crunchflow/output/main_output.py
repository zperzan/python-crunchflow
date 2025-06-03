"""
Module for loading and parsing the main CrunchFlow output files.

Defines classes to read and organize information within the main CrunchFlow output file.
If your main input file is run_name.in, this file is named run_name.out. As of now, only
the speciation blocks of each geochemical condition are read in, but future versions may
expand the information read from this file.
"""

from pathlib import Path

from .speciation import Speciation


class SpeciationBlockCollection(dict):
    """
    A specialized dictionary to hold multiple Speciation objects keyed by condition name.

    Provides a custom string representation listing all contained geochemical conditions.
    """

    def __str__(self):
        """Return a string listing all geochemical conditions in the collection."""
        print_str = f"SpeciationBlockCollection with {len(self)} conditions: {', '.join(self.keys())}"
        return print_str


class MainOutputFile:
    """
    Parses a CrunchFlow main output file to extract speciation blocks.

    Attributes
    ----------
    filepath : Path
        The path to the CrunchFlow output file.
    speciation_blocks : SpeciationBlockCollection
        Dictionary-like collection of Speciation objects indexed by condition name.

    Methods
    -------
    _load_file()
        Reads the file and populates the speciation_blocks attribute.
    """

    def __init__(self, filepath):
        """
        Initialize the MainOutputFile by reading the specified file.

        Parameters
        ----------
        filepath : str or Path
            The path to the CrunchFlow output file.
        """
        self.filepath = Path(filepath)
        self.speciation_blocks = SpeciationBlockCollection()

        self._load_file()

    def _load_file(self):
        """
        Parse the CrunchFlow output file and extract speciation blocks.

        Detects when a geochemical condition starts and ends, and uses that to construct
        Speciation objects that are stored in speciation_blocks.
        """
        if not self.filepath.exists():
            raise FileNotFoundError(f"No such file: {self.filepath}")

        block_lines = []
        inside_speciation_section = False
        inside_block = False

        with open(self.filepath, "r", encoding="utf8", errors="ignore") as f:
            for line in f:
                # Determine whether we're in the SPECIATION section
                if "SPECIATION OF GEOCHEMICAL CONDITIONS" in line:
                    inside_speciation_section = True
                elif "SUCCESSFULLY COMPLETED" in line or "DATABASE SWEEP SPECIFIED" in line:
                    inside_speciation_section = False

                    # Append final block if it exists
                    if block_lines:
                        block_text = "\n".join(block_lines)
                        if len(block_text.strip()) != 0:
                            self.speciation_blocks[condition_name] = Speciation(block_text)

                if inside_speciation_section:
                    # Use GEOCHEMICAL CONDITION to identify start of a block
                    if "GEOCHEMICAL CONDITION:" in line:
                        condition_name = line.split(":")[-1].strip()
                        inside_block = True
                        block_lines = []
                        continue
                    # Use asterisks to identify end of previous block
                    elif "*************************" in line:
                        inside_block = False
                        block_text = "\n".join(block_lines)
                        if len(block_text.strip()) != 0:
                            self.speciation_blocks[condition_name] = Speciation(block_text)

                    if inside_block and line.strip() != 0:
                        block_lines.append(line)
