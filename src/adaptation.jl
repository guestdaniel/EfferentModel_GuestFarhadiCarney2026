export adapt_pla, adapt_ea, adapt_ea_iir, adapt_ea_iir_parallel

function adapt_pla(x, α, β; fs=10e3)
    y = zeros(size(x))
    I = zeros(size(x))
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - α * I[n-1])
        end

        # Compute I[n]
        for j in 1:n
            I[n] += y[j] * 1/fs / ((n-j)*1/fs + β)
        end
    end
    return y, I
end

function adapt_ea_iir_parallel(x, α, τs; fs=10e3)
    y = zeros(length(x))
    I = zeros(length(x), length(τs))
    I_comb = zeros(length(x))
    ds = @. exp(-(1/fs)/τs)
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - α * I_comb[n-1])
        end

        # Compute I[n]
        for j in eachindex(ds)
            if n == 1
                I[n, j] = (1-ds[j]) * 1e-1 * y[n]
            else
                I[n, j] = (1-ds[j]) * 1e-1 * y[n] + ds[j] * I[n-1, j]
            #    I[n, j] = (1-ds[j]) * 1e-1 * y[n] + ds[j] * I_comb[n-1]
            end
        end
        I_comb[n] = sum(I[n, :])
    end
    return y, I_comb, I
end

function adapt_ea(x, τₐ, τₑ; fs=10e3)
    y = zeros(size(x))
    I = zeros(size(x))
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - (1/τₐ) * I[n-1])
        end

        # Compute I[n]
        for j in 1:n
            I[n] += y[j] * 1/fs * exp( ((j-n) * (1/fs)) * (1/τₑ))
        end
    end
    return y, I
end

function adapt_ea_iir(x, τₐ, τ; fs=10e3)
    y = zeros(size(x))
    I = zeros(size(x))
    d = exp(-(1/fs)/τ)
    for n in eachindex(x)
        # Apply PLA
        if n == 1
            y[n] = max(0.0, x[n])
        else
            y[n] = max(0.0, x[n] - (1/τₐ) * I[n-1])
        end

        # Compute I[n]
        if n == 1
            I[n] = (1-d) * 1e-1 * y[n]
        else
            I[n] = (1-d) * 1e-1 * y[n] + d * I[n-1]
        end
    end
    return y, I
end

