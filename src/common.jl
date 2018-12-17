const spectroscopic = "spdfghiklmnoqrtuv"
spectroscopic_label(ℓ) =
    ℓ + 1 ≤ length(spectroscopic) ? spectroscopic[ℓ+1] : "[$(ℓ)]"

# Nicer string representation for rationals
rs(r::Number) = "$(r)"
rs(r::Rational) = "$(numerator(r))/$(denominator(r))"
