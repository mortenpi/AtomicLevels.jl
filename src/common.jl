spectroscopic = "spdfghiklmnoqrtuv"
spectroscopic_label(ℓ) =
    ℓ + 1 ≤ length(spectroscopic) ? spectroscopic[ℓ+1] : "[$(ℓ)]"
