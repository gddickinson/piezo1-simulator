"""The 35 harvested candidates, read by hand, one verdict each.

Round 45 extracted 35 substitutions from the open-access corpus that are not in
the curated set, each with the sentence it came from. Round 57's question was
precise: **how many have a direction recoverable by a human from the sentence
alone?** That number decides whether another test is ever possible.

Every sentence was read. The answer is **five**, and all five say the same
thing — "non-functional" — which is loss of *channel* function rather than the
loss-of-function-in-disease the curated set records.

The rest fail for reasons worth separating, because they are different problems:

===================  ==  ====================================================
category              n  what it is
===================  ==  ====================================================
loss_of_function      5  the sentence states the channel does not work
chemical_only         4  a direction for agonist response, not for mechanics
conductance_only      5  a conductance change, which is not a direction
no_phenotype         17  a construct list, a figure legend, or a location
sequence_variant      3  a clone's difference from the reference, not a mutant
wrong_protein         1  the sentence is about STOML3, not PIEZO1
===================  ==  ====================================================

**The one that matters most is `wrong_protein`.** V190P is a STOML3 mutation
from a paper that studies STOML3 and PIEZO1 together. The harvest's wild-type
gate passed it because position 190 is valine in PIEZO1 too — so the gate,
which rejects 23% of raw hits and is the reason to trust the rest, cannot catch
a substitution that is real but belongs to another protein. That is a class of
error no residue-identity check can see, and it is recorded rather than quietly
filtered.

**And `conductance_only` is the same question Round 54 left open.** A change in
unitary conductance is a measured functional effect that is not a gain or loss
of mechanosensitivity. Whether one may stand for the other is a scientific
decision, not a curation one, so these are held here rather than admitted.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["Verdict", "CURATION", "CATEGORIES", "by_category", "directional",
           "summary"]

#: Why a candidate does or does not yield a direction. Kept explicit because
#: "35 candidates" and "5 usable" are very different numbers to plan with.
CATEGORIES = {
    "loss_of_function": "the sentence states the channel is non-functional",
    "chemical_only": "a direction for agonist response, not for mechanics",
    "conductance_only": "a conductance change, which is not a direction",
    "no_phenotype": "a construct list, figure legend or location statement",
    "sequence_variant": "a clone's difference from the reference sequence",
    "wrong_protein": "the sentence describes a different protein",
}


@dataclass(frozen=True)
class Verdict:
    """One candidate, read by hand."""

    label: str                # human numbering
    category: str
    direction: str | None     # "LoF" only where the sentence supports it
    basis: str                # the phrase the verdict rests on

    def __post_init__(self) -> None:
        if self.category not in CATEGORIES:
            raise ValueError(f"unknown category {self.category!r}")
        if self.direction is not None and self.category != "loss_of_function":
            raise ValueError(
                f"{self.label}: a direction may only come from a stated "
                f"functional outcome, not from {self.category}")


CURATION: tuple = (
    # --- the five that carry a direction -------------------------------
    Verdict("D1959A", "loss_of_function", "LoF",
            "D1975A and D2034A were non-functional mutants"),
    Verdict("D2018A", "loss_of_function", "LoF",
            "D1975A and D2034A were non-functional mutants"),
    Verdict("L2131A", "loss_of_function", "LoF",
            "resulted either in non-functional channels (P2129A, L2131A, "
            "R2135A and W2140A)"),
    Verdict("R2119A", "loss_of_function", "LoF",
            "non-functional channels (... R2135A ...)"),
    Verdict("W2124A", "loss_of_function", "LoF",
            "non-functional channels (... W2140A)"),

    # --- agonist response only, and two of them are double mutants -----
    Verdict("L1347A", "chemical_only", None,
            "the L1342A/L1345A DOUBLE mutant prevents chemical activation and "
            "impairs poking, but maintains stretch sensitivity — the single "
            "mutant's own phenotype is not stated"),
    Verdict("L1350A", "chemical_only", None,
            "same double-mutant sentence; the single mutant is not separated"),
    Verdict("A1718G", "chemical_only", None,
            "retains partial Yoda1 sensitivity and remains activatable by "
            "Jedi2 — a chemical-agonist phenotype"),
    Verdict("A1718V", "chemical_only", None,
            "retains partial Yoda1 sensitivity — chemical, not mechanical"),

    # --- conductance changes: measured, but not a direction ------------
    Verdict("F2114A", "conductance_only", None,
            "unitary conductance similar to WT — an explicit NO-change result"),
    Verdict("D2123A", "conductance_only", None,
            "slightly increased single-channel conductance in divalent-free "
            "conditions"),
    Verdict("D2128A", "conductance_only", None,
            "slightly increased single-channel conductance"),
    Verdict("V2116A", "conductance_only", None,
            "the sentence is a fragment of a conductance list ('7 pS, "
            "V2132A; 59.') and the parsed measurement is a truncation artefact"),
    Verdict("L2118A", "conductance_only", None,
            "same truncated conductance list; the value is not recoverable"),

    # --- no phenotype in the sentence ----------------------------------
    Verdict("L1347D", "no_phenotype", None,
            "a list of constructs made, with no outcome"),
    Verdict("L1347S", "no_phenotype", None, "'we mutated L1342 and L1345 to different amino acids, resulting in "
            "the following mutants' — constructs made, no outcome"),
    Verdict("L1350D", "no_phenotype", None, "'we mutated L1342 and L1345 to different amino acids, resulting in "
            "the following mutants' — constructs made, no outcome"),
    Verdict("L1350S", "no_phenotype", None, "'we mutated L1342 and L1345 to different amino acids, resulting in "
            "the following mutants' — constructs made, no outcome"),
    Verdict("Q1349A", "no_phenotype", None,
            "did not affect Jedi1- or Yoda1-induced responses — a null result, "
            "which is informative and is not a direction"),
    Verdict("A1718I", "no_phenotype", None,
            "from a figure legend listing plot symbols"),
    Verdict("A1718L", "no_phenotype", None, "from a figure legend listing plot symbols for each tested mutant"),
    Verdict("A2078D", "no_phenotype", None, "from a figure legend listing plot symbols for each tested mutant"),
    Verdict("A2078F", "no_phenotype", None, "from a figure legend listing plot symbols for each tested mutant"),
    Verdict("A2078V", "no_phenotype", None, "from a figure legend listing plot symbols for each tested mutant"),
    Verdict("D1959N", "no_phenotype", None,
            "'we used D1975N and D2034N for further characterization' — the "
            "outcome is in an adjacent sentence, not this one"),
    Verdict("D2018N", "no_phenotype", None,
            "'we used D1975N and D2034N for further characterization' — the "
            "outcome is in an adjacent sentence, not this one"),
    Verdict("R1926Q", "no_phenotype", None,
            "mutations ... are scattered throughout the channel — locational"),
    Verdict("R2285H", "no_phenotype", None,
            "the same 'scattered throughout the channel' sentence — locational"),
    Verdict("R2119A", "no_phenotype", None,
            "the botellosmith2019 hit for this residue is a figure legend; the "
            "coste2015 hit for the same residue is the one that carries LoF"),
    Verdict("L1347A", "no_phenotype", None,
            "the second source for this residue adds no separate phenotype"),
    Verdict("L1350A", "no_phenotype", None,
            "the second source for this residue repeats the double-mutant "
            "sentence and adds no separate phenotype"),

    # --- clone sequence differences, not tested mutants -----------------
    Verdict("V250A", "sequence_variant", None,
            "this clone differs from the NCBI sequence at ... — a sequencing "
            "note about the construct, not a mutant that was tested"),
    Verdict("V394L", "sequence_variant", None,
            "the same clone-versus-NCBI sequencing note, not a tested mutant"),
    Verdict("R407G", "sequence_variant", None,
            "the same clone-versus-NCBI sequencing note, not a tested mutant"),

    # --- not PIEZO1 at all ---------------------------------------------
    Verdict("V190P", "wrong_protein", None,
            "'we introduced a mutation into STOML3 at the orthologous position "
            "(V190P)' — STOML3, not PIEZO1. The wild-type gate passed it "
            "because position 190 is valine in PIEZO1 as well."),
)


def by_category(category: str) -> list:
    return [v for v in CURATION if v.category == category]


def directional() -> list:
    """Candidates a human can assign a direction to from the sentence alone."""
    return [v for v in CURATION if v.direction is not None]


def summary() -> dict:
    counts = {name: len(by_category(name)) for name in CATEGORIES}
    return {"total": len(CURATION), "by_category": counts,
            "directional": len(directional()),
            "distinct_directional": len({v.label for v in directional()})}
