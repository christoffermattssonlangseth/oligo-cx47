from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class QuestionSpec:
    identifier: str
    short_title: str
    notebook_slug: str
    full_question: str
    rationale: str
    panels: tuple[str, ...]
    recommended_obs_keywords: tuple[str, ...]
    dataset_ids: tuple[str, ...]
    optional: bool = False


QUESTION_BANK: dict[str, QuestionSpec] = {
    "q1": QuestionSpec(
        identifier="q1",
        short_title="Panglial Network Remodeling",
        notebook_slug="01_panglial_network_remodeling",
        full_question="How does the panglial network spatiotemporally remodel during demyelination and repair?",
        rationale=(
            "Track oligodendrocyte and astrocyte connexins across lesion and repair states, "
            "with emphasis on Cx47, Cx32, and Cx43 re-establishment during remyelination."
        ),
        panels=("panglial_connexins", "oligodendrocyte_identity", "astrocyte_identity", "ion_homeostasis"),
        recommended_obs_keywords=("time", "day", "week", "condition", "cell", "annotation", "lesion", "region"),
        dataset_ids=("lpc_remyelination", "ms_heterogeneity", "eae_ms_adaptive_oligo"),
    ),
    "q2": QuestionSpec(
        identifier="q2",
        short_title="Cx47 and Mitochondrial Recovery",
        notebook_slug="02_cx47_mitochondrial_recovery",
        full_question="May Cx47 mark or coordinate the mitochondrial metabolic transition of oligodendrocytes toward remyelination?",
        rationale=(
            "Relate Cx47 re-expression to OXPHOS, mitochondrial biogenesis, transport, and "
            "metabolic coupling signatures across LPC remyelination time points."
        ),
        panels=(
            "panglial_connexins",
            "mitochondrial_oxphos",
            "mitochondrial_membrane_transport",
            "mitochondrial_dynamics",
            "mitochondrial_biogenesis",
            "metabolic_coupling",
        ),
        recommended_obs_keywords=("time", "day", "week", "condition", "cell", "cluster", "oligo"),
        dataset_ids=("lpc_remyelination", "ms_heterogeneity"),
    ),
    "q3": QuestionSpec(
        identifier="q3",
        short_title="ER Stress and Reactive Glia",
        notebook_slug="03_er_stress_and_reactive_glia",
        full_question="How do ER stress and reactive glial states influence panglial network remodeling?",
        rationale=(
            "Test whether ER stress, inflammatory cytokines, reactive astrocytes, and activated "
            "microglia are enriched early after demyelination and decline during Cx47 recovery."
        ),
        panels=(
            "panglial_connexins",
            "er_stress_upr",
            "inflammatory_cytokines",
            "microglial_activation",
            "astrocyte_identity",
        ),
        recommended_obs_keywords=("time", "condition", "cell", "annotation", "microglia", "astro", "cluster"),
        dataset_ids=("lpc_remyelination", "eae_ms_adaptive_oligo", "cx47_cko_eae"),
    ),
    "q4": QuestionSpec(
        identifier="q4",
        short_title="Transglial Mitochondrial Communication",
        notebook_slug="04_transglial_mitochondrial_communication",
        full_question="Does Cx47 regulate transglial mitochondrial communication during neuroinflammation?",
        rationale=(
            "Explore metabolic coupling, mitochondrial trafficking, and stress programs across "
            "neighboring glial populations as a proxy for transglial mitochondrial communication."
        ),
        panels=(
            "panglial_connexins",
            "mitochondrial_trafficking",
            "mitophagy_stress",
            "metabolic_coupling",
            "microglial_activation",
            "astrocyte_identity",
        ),
        recommended_obs_keywords=("time", "condition", "cell", "neighbor", "region", "sample"),
        dataset_ids=("lpc_remyelination", "cx47_cko_eae", "cuprizone_trem2ko"),
        optional=True,
    ),
    "q5": QuestionSpec(
        identifier="q5",
        short_title="Lipid Metabolism and Glymphatic Links",
        notebook_slug="05_lipid_metabolism_and_glymphatic_links",
        full_question="Does Cx47 influence lipid metabolism and clearance pathways linking panglia to meningeal glymphatic systems during neuroinflammation?",
        rationale=(
            "Probe whether lipid handling, lymphatic endothelial, and barrier-associated pathways "
            "co-vary with connexin state during inflammatory demyelination."
        ),
        panels=(
            "panglial_connexins",
            "lipid_metabolism",
            "lymphatic_endothelial",
            "lymphatic_barrier",
            "bbb_related",
        ),
        recommended_obs_keywords=("time", "condition", "cell", "endothelial", "lymph", "region", "sample"),
        dataset_ids=("lpc_remyelination", "ms_heterogeneity", "cuprizone_trem2ko"),
        optional=True,
    ),
}

ORDERED_QUESTION_IDS: tuple[str, ...] = tuple(QUESTION_BANK)
