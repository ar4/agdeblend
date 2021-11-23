module FakeMPI
  Comm = Nothing
  MPI_Comm = nothing
end

module AGDeblend

haveMPI = try
  using MPI
  true
catch
  false
end
  
if (haveMPI)
  using MPI
else
  const MPI = FakeMPI
end

AGDLive = Cint(0)
AGDBad = Cint(1)
AGDMissing = Cint(2)
AGDBlendSum = Cint(0)
AGDBlendMean = Cint(1)
AGDBlendOverwrite = Cint(2)

export blend, deblend, AGDLive, AGDBad, AGDMissing, AGDBlendSum,
    AGDBlendMean, AGDBlendOverwrite

function deblend(volumes::Array{<:Integer},
                 window_shapes::Array{<:Array{<:Integer}},
                 coords::Array{<:Array{<:Integer}},
                 shottimes::Array{<:Array{<:Integer}},
                 channels::Array{<:Array{<:Integer}},
                 trace_types::Array{<:Array{<:Integer}},
                 initial_factor::AbstractFloat,
                 n_its::Integer, print_freq::Integer,
                 data::Array{<:Array{<:AbstractFloat}};
                 wavelet_idxs::Union{Array{<:Array{<:Integer}},
                                     Nothing}=nothing,
                 wavelets::Union{Array{<:Array{<:AbstractFloat}},
                                 Nothing}=nothing,
                 comm::Union{MPI.Comm, Nothing}=nothing)::Nothing
    dtype = eltype(eltype(data))
    if !(dtype in [Cfloat, Cdouble])
        error("data must be equivalent to Cfloat or Cdouble")
    end
    n_patches = size(data, 1)
    shapes = [Array{Cint}(collect(size(patch_data)))
              for patch_data in data]
    if wavelets !== nothing
        if wavelet_idxs === nothing
            error("wavelets and wavelet_idxs must either both be nothing or "
                  * "neither")
        end
        wavelet_lengths = [Cint(length(wavelet)) for wavelet in wavelets]
        wavelet_idxs = Array{Array{Cint}}(wavelet_idxs)
        wavelets = Array{Array{dtype}}(wavelets)
    else
        if wavelet_idxs !== nothing
            error("wavelets and wavelet_idxs must either both be nothing or "
                  * "neither")
        end
        wavelet_lengths = C_NULL
        wavelet_idxs = C_NULL
        wavelets = C_NULL
    end
    
    # Check that the number of patches is the same for all variables
    if length(volumes) != n_patches
        error("volumes must have the same length as data")
    end
    check_n_patches(n_patches, coords, "coords")
    check_n_patches(n_patches, shottimes, "shottimes")
    check_n_patches(n_patches, channels, "channels")
    check_n_patches(n_patches, trace_types, "trace_types")
    if wavelet_idxs != C_NULL
        check_n_patches(n_patches, wavelet_idxs, "wavelet_idxs")
    end
    
    # Check that the number of traces in each patch is the same for
    # all variables
    n_traces = [prod(size(patch_data)[2:end]) for patch_data in data]
    check_n_traces(n_traces, shottimes, "shottimes")
    check_n_traces(n_traces, channels, "channels")
    check_n_traces(n_traces, trace_types, "trace_types")
    if wavelet_idxs != C_NULL
        check_n_traces(n_traces, wavelet_idxs, "wavelet_idxs")
    end
    
    # Check wavelet_idxs match wavelets
    if wavelet_idxs != C_NULL
        check_wavelets(wavelet_idxs, wavelets, trace_types)
    end

    # Extract the number of dimensions in each volume
    n_volumes = maximum(volumes) + 1
    n_dims = zeros(Cint, n_volumes)
    for (volume_idx, patch_data) in zip(volumes, data)
      n_dims[volume_idx + 1] = Cint(ndims(patch_data))
    end

    # Flip arrays to switch to row-major
    window_shapes = [window_shape[end:-1:1] for window_shape in window_shapes]
    coords = [coord[end:-1:1] for coord in coords]
    shapes = [patch_shape[end:-1:1] for patch_shape in shapes]
    
    ret = deblend_c(n_patches, Array{Cint}(volumes), Array{Cint}(n_dims),
                    Array{Array{Cint}}(window_shapes),
                    Array{Array{Cint}}(coords),
                    Array{Array{Cint}}(shapes),
                    Array{Array{Clong}}(shottimes),
                    Array{Array{Cint}}(channels),
                    Array{Array{Cint}}(trace_types),
                    wavelet_lengths, wavelet_idxs, wavelets, initial_factor,
                    n_its, print_freq, comm, data)
    if ret != 0
        error("ERROR in deblend")
    end

    return nothing
end

# Single precision, no MPI
function deblend_c(n_patches, volumes::Array{Cint}, n_dims::Array{Cint},
                   window_shapes::Array{<:Array{Cint}},
                 coords::Array{<:Array{Cint}},
                 shapes::Array{<:Array{Cint}},
                 shottimes::Array{<:Array{Clong}},
                 channels::Array{<:Array{Cint}},
                 trace_types::Array{<:Array{Cint}},
                 wavelet_lengths::Union{Array{Cint}, Ptr{Nothing}},
                 wavelet_idxs::Union{Array{<:Array{Cint}}, Ptr{Nothing}},
                 wavelets::Union{Array{<:Array{Cfloat}}, Ptr{Nothing}},
                 initial_factor, n_its, print_freq, comm::Nothing,
                 data::Array{<:Array{Cfloat}})
    return ccall((:agd_deblend, "libagdeblend"), Cint,
                 (Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Cint}}, Ref{Ptr{Cint}},
                  Ref{Ptr{Cint}}, Ref{Ptr{Clong}}, Ref{Ptr{Cint}},
                  Ref{Ptr{Cint}}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cfloat}},
                  Cfloat, Cint, Cint, Ref{Ptr{Cfloat}}),
                 n_patches, volumes, n_dims,
                 window_shapes, coords,
                 shapes, shottimes,
                 channels, trace_types,
                 wavelet_lengths, wavelet_idxs, wavelets, initial_factor, n_its,
                 print_freq, data)
end

# Double precision, no MPI
function deblend_c(n_patches, volumes::Array{Cint}, n_dims::Array{Cint},
                   window_shapes::Array{<:Array{Cint}},
                 coords::Array{<:Array{Cint}},
                 shapes::Array{<:Array{Cint}},
                 shottimes::Array{<:Array{Clong}},
                 channels::Array{<:Array{Cint}},
                 trace_types::Array{<:Array{Cint}},
                 wavelet_lengths::Union{Array{Cint}, Ptr{Nothing}},
                 wavelet_idxs::Union{Array{<:Array{Cint}}, Ptr{Nothing}},
                 wavelets::Union{Array{<:Array{Cdouble}}, Ptr{Nothing}},
                 initial_factor, n_its, print_freq, comm::Nothing,
                 data::Array{<:Array{Cdouble}})
    return ccall((:agd_deblend, "libagdeblend_double"), Cint,
                 (Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Cint}}, Ref{Ptr{Cint}},
                  Ref{Ptr{Cint}}, Ref{Ptr{Clong}}, Ref{Ptr{Cint}},
                  Ref{Ptr{Cint}}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}},
                  Cdouble, Cint, Cint, Ref{Ptr{Cdouble}}),
                 n_patches, volumes, n_dims,
                 window_shapes, coords,
                 shapes, shottimes,
                 channels, trace_types,
                 wavelet_lengths, wavelet_idxs, wavelets, initial_factor, n_its,
                 print_freq, data)
end

if haveMPI
# Single precision, MPI
  function deblend_c(n_patches, volumes::Array{Cint}, n_dims::Array{Cint},
                     window_shapes::Array{<:Array{Cint}},
                   coords::Array{<:Array{Cint}},
                   shapes::Array{<:Array{Cint}},
                   shottimes::Array{<:Array{Clong}},
                   channels::Array{<:Array{Cint}},
                   trace_types::Array{<:Array{Cint}},
                   wavelet_lengths::Union{Array{Cint}, Ptr{Nothing}},
                   wavelet_idxs::Union{Array{<:Array{Cint}}, Ptr{Nothing}},
                   wavelets::Union{Array{<:Array{Cfloat}}, Ptr{Nothing}},
                   initial_factor, n_its, print_freq, comm::MPI.Comm,
                   data::Array{<:Array{Cfloat}})
      return ccall((:agd_deblend, "libagdeblend_mpi"), Cint,
                   (Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Cint}}, Ref{Ptr{Cint}},
                    Ref{Ptr{Cint}}, Ref{Ptr{Clong}}, Ref{Ptr{Cint}},
                    Ref{Ptr{Cint}}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cfloat}},
                    Cfloat, Cint, Cint, MPI.MPI_Comm, Ref{Ptr{Cfloat}}),
                   n_patches, volumes, n_dims,
                   window_shapes, coords,
                   shapes, shottimes,
                   channels, trace_types,
                   wavelet_lengths, wavelet_idxs, wavelets, initial_factor, n_its,
                   print_freq, comm, data)
  end
  
  # Double precision, MPI
  function deblend_c(n_patches, volumes::Array{Cint}, n_dims::Array{Cint},
                     window_shapes::Array{<:Array{Cint}},
                   coords::Array{<:Array{Cint}},
                   shapes::Array{<:Array{Cint}},
                   shottimes::Array{<:Array{Clong}},
                   channels::Array{<:Array{Cint}},
                   trace_types::Array{<:Array{Cint}},
                   wavelet_lengths::Union{Array{Cint}, Ptr{Nothing}},
                   wavelet_idxs::Union{Array{<:Array{Cint}}, Ptr{Nothing}},
                   wavelets::Union{Array{<:Array{Cdouble}}, Ptr{Nothing}},
                   initial_factor, n_its, print_freq, comm::MPI.Comm,
                   data::Array{<:Array{Cdouble}})
      return ccall((:agd_deblend, "libagdeblend_double_mpi"), Cint,
                   (Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Cint}}, Ref{Ptr{Cint}},
                    Ref{Ptr{Cint}}, Ref{Ptr{Clong}}, Ref{Ptr{Cint}},
                    Ref{Ptr{Cint}}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}},
                    Cdouble, Cint, Cint, MPI.MPI_Comm, Ref{Ptr{Cdouble}}),
                   n_patches, volumes, n_dims,
                   window_shapes, coords,
                   shapes, shottimes,
                   channels, trace_types,
                   wavelet_lengths, wavelet_idxs, wavelets, initial_factor, n_its,
                   print_freq, comm, data)
  end
end

function blend(shottimes::Array{<:Array{<:Integer}},
               channels::Array{<:Array{<:Integer}},
               trace_types::Array{<:Array{<:Integer}},
               data::Array{<:Array{<:AbstractFloat}},
               blend_mode::Integer;
               taper_length::Integer=0,
               n_times_out::Union{Array{<:Integer}, Nothing}=nothing,
               shottimes_out::Union{Array{<:Array{<:Integer}}, Nothing}=nothing,
               channels_out::Union{Array{<:Array{<:Integer}}, Nothing}=nothing,
               trace_types_out::Union{Array{<:Array{<:Integer}},
                                      Nothing}=nothing,
               comm::Union{MPI.Comm, Nothing}=nothing,
               data_out::Union{Array{<:Array{<:AbstractFloat}},
                               Nothing}=nothing)::Array{<:Array{<:AbstractFloat}}
    dtype = eltype(eltype(data))
    if !(dtype in [Cfloat, Cdouble])
        error("data must be equivalent to Cfloat or Cdouble")
    end
    if !(blend_mode in [AGDBlendSum, AGDBlendMean, AGDBlendOverwrite])
        error("Invalid blend_mode")
    end
    shottimes_out = shottimes_out === nothing ? shottimes : shottimes_out
    channels_out = channels_out === nothing ? channels : channels_out
    trace_types_out = trace_types_out === nothing ? trace_types : trace_types_out
    
    n_patches = size(data, 1)
    n_traces = [Cint(prod(size(patch_data)[2:end])) for patch_data in data]
    n_times = [Cint(size(patch_data)[1]) for patch_data in data]
    
    # If data_out is not provided, use shottimes_out and n_times_out to get
    # the output size
    if data_out !== nothing
        if length(data_out) < 1
            error("data_out must have at least one patch")
        end
        if eltype(eltype(data_out)) != dtype
            error("data_out must have the same data type as data")
        end
        n_patches_out = size(data_out, 1)
        n_traces_out = [Cint(prod(size(patch_data)[2:end]))
                        for patch_data in data_out]
        if n_times_out === nothing
            n_times_out = [size(patch_data)[end] for patch_data in data_out]
        else
            if n_times_out != [size(patch_data)[end] for patch_data in data_out]
                error("data_out does not match n_times_out")
            end
        end
        n_times_out = Array{Cint}(n_times_out)
    else # data_out not provided, so use shottimes_out and n_times_out/n_times
        n_patches_out = size(shottimes_out, 1)
        n_traces_out = [Cint(prod(size(patch_shottimes)))
                        for patch_shottimes in shottimes_out]
        if n_times_out === nothing
            n_times_out = n_times
        else
            n_times_out = Array{Cint}(n_times_out)
        end
        data_out = [zeros(dtype, append!([patch_n_times],
                                         collect(size(patch_shottimes)))...)
                    for (patch_shottimes, patch_n_times) in
                    zip(shottimes_out, n_times_out)]
    end
    
    if length(n_times_out) != n_patches_out
        error("n_times_out is not of the expected length")
    end
    
    # Check that the number of patches is the same for all variables
    check_n_patches(n_patches, shottimes, "shottimes")
    check_n_patches(n_patches, channels, "channels")
    check_n_patches(n_patches, trace_types, "trace_types")
    check_n_patches(n_patches_out, shottimes_out, "shottimes_out")
    check_n_patches(n_patches_out, channels_out, "channels_out")
    check_n_patches(n_patches_out, trace_types_out, "trace_types_out")
    
    # Check that the number of traces in each patch is the same for
    # all variables
    check_n_traces(n_traces, shottimes, "shottimes")
    check_n_traces(n_traces, channels, "channels")
    check_n_traces(n_traces, trace_types, "trace_types")
    check_n_traces(n_traces_out, shottimes_out, "shottimes_out")
    check_n_traces(n_traces_out, channels_out, "channels_out")
    check_n_traces(n_traces_out, trace_types_out, "trace_types_out")
    
    ret = blend_c(n_patches, n_traces, n_times,
                    Array{Array{Clong}}(shottimes),
                    Array{Array{Cint}}(channels), Array{Array{Cint}}(trace_types),
                    data, blend_mode, taper_length, n_patches_out,
                    n_traces_out, n_times_out, Array{Array{Clong}}(shottimes_out),
                    Array{Array{Cint}}(channels_out),
                    Array{Array{Cint}}(trace_types_out), comm, data_out)
    if ret != 0
        error("ERROR in blend")
    end

    return data_out
end

# Single precision, no MPI
function blend_c(n_patches, n_traces::Array{Cint}, n_times::Array{Cint},
                 shottimes::Array{<:Array{Clong}},
                 channels::Array{<:Array{Cint}},
                 trace_types::Array{<:Array{Cint}},
                 data::Array{<:Array{Cfloat}},
                 blend_mode, taper_length, n_patches_out,
                 n_traces_out::Array{Cint}, n_times_out::Array{Cint},
                 shottimes_out::Array{<:Array{Clong}},
                 channels_out::Array{<:Array{Cint}},
                 trace_types_out::Array{<:Array{Cint}},
                 comm::Nothing, data_out::Array{<:Array{Cfloat}})
    return ccall((:agd_blend, "libagdeblend"), Cint,
                 (Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Clong}}, Ref{Ptr{Cint}},
                  Ref{Ptr{Cint}}, Ref{Ptr{Cfloat}}, Cint,
                  Cint, Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Clong}},
                  Ref{Ptr{Cint}}, Ref{Ptr{Cint}}, Ref{Ptr{Cfloat}}),
                 n_patches, n_traces, n_times, shottimes, channels, trace_types,
                 data, blend_mode, taper_length, n_patches_out, n_traces_out,
                 n_times_out, shottimes_out, channels_out, trace_types_out,
                 data_out)
end

# Double precision, no MPI
function blend_c(n_patches, n_traces::Array{Cint}, n_times::Array{Cint},
                 shottimes::Array{<:Array{Clong}},
                 channels::Array{<:Array{Cint}},
                 trace_types::Array{<:Array{Cint}},
                 data::Array{<:Array{Cdouble}},
                 blend_mode, taper_length, n_patches_out,
                 n_traces_out::Array{Cint}, n_times_out::Array{Cint},
                 shottimes_out::Array{<:Array{Clong}},
                 channels_out::Array{<:Array{Cint}},
                 trace_types_out::Array{<:Array{Cint}},
                 comm::Nothing, data_out::Array{<:Array{Cdouble}})
    return ccall((:agd_blend, "libagdeblend_double"), Cint,
                 (Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Clong}}, Ref{Ptr{Cint}},
                  Ref{Ptr{Cint}}, Ref{Ptr{Cdouble}}, Cint,
                  Cint, Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Clong}},
                  Ref{Ptr{Cint}}, Ref{Ptr{Cint}}, Ref{Ptr{Cdouble}}),
                 n_patches, n_traces, n_times, shottimes, channels, trace_types,
                 data, blend_mode, taper_length, n_patches_out, n_traces_out,
                 n_times_out, shottimes_out, channels_out, trace_types_out,
                 data_out)
end

if haveMPI
  # Single precision, MPI
  function blend_c(n_patches, n_traces::Array{Cint}, n_times::Array{Cint},
                   shottimes::Array{<:Array{Clong}},
                   channels::Array{<:Array{Cint}},
                   trace_types::Array{<:Array{Cint}},
                   data::Array{<:Array{Cfloat}},
                   blend_mode, taper_length, n_patches_out,
                   n_traces_out::Array{Cint}, n_times_out::Array{Cint},
                   shottimes_out::Array{<:Array{Clong}},
                   channels_out::Array{<:Array{Cint}},
                   trace_types_out::Array{<:Array{Cint}},
                   comm::MPI.Comm, data_out::Array{<:Array{Cfloat}})
      return ccall((:agd_blend, "libagdeblend_mpi"), Cint,
                   (Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Clong}}, Ref{Ptr{Cint}},
                    Ref{Ptr{Cint}}, Ref{Ptr{Cfloat}}, Cint,
                    Cint, Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Clong}},
                    Ref{Ptr{Cint}}, Ref{Ptr{Cint}}, MPI.MPI_Comm,
                    Ref{Ptr{Cfloat}}),
                   n_patches, n_traces, n_times, shottimes, channels, trace_types,
                   data, blend_mode, taper_length, n_patches_out, n_traces_out,
                   n_times_out, shottimes_out, channels_out, trace_types_out,
                   comm, data_out)
  end

  # Double precision, MPI
  function blend_c(n_patches, n_traces::Array{Cint}, n_times::Array{Cint},
                   shottimes::Array{<:Array{Clong}},
                   channels::Array{<:Array{Cint}},
                   trace_types::Array{<:Array{Cint}},
                   data::Array{<:Array{Cdouble}},
                   blend_mode, taper_length, n_patches_out,
                   n_traces_out::Array{Cint}, n_times_out::Array{Cint},
                   shottimes_out::Array{<:Array{Clong}},
                   channels_out::Array{<:Array{Cint}},
                   trace_types_out::Array{<:Array{Cint}},
                   comm::MPI.Comm, data_out::Array{<:Array{Cdouble}})
      return ccall((:agd_blend, "libagdeblend_double_mpi"), Cint,
                   (Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Clong}}, Ref{Ptr{Cint}},
                    Ref{Ptr{Cint}}, Ref{Ptr{Cdouble}}, Cint,
                    Cint, Cint, Ref{Cint}, Ref{Cint}, Ref{Ptr{Clong}},
                    Ref{Ptr{Cint}}, Ref{Ptr{Cint}}, MPI.MPI_Comm,
                    Ref{Ptr{Cdouble}}),
                   n_patches, n_traces, n_times, shottimes, channels, trace_types,
                   data, blend_mode, taper_length, n_patches_out, n_traces_out,
                   n_times_out, shottimes_out, channels_out, trace_types_out,
                   comm, data_out)
  end
end

function check_n_patches(n_patches, var, name)
    if length(var) != n_patches
        error(name * " must have the same length as data")
    end
end

function get_n_traces(var)
    return [prod(size(patch_var)) for patch_var in var]
end

function check_n_traces(n_traces, var, name)
    if get_n_traces(var) != n_traces
        error(name * " has an incorrect number of traces")
    end
end

function check_wavelets(wavelet_idxs, wavelets, trace_types)
    min_wavelet_idx = 0
    max_wavelet_idx = 0
    for (patch_wavelet_idxs, patch_trace_types) in zip(wavelet_idxs,
                                                       trace_types)
        flat_patch_wavelet_idxs = collect(Iterators.flatten(patch_wavelet_idxs))
        flat_patch_trace_types = collect(Iterators.flatten(patch_trace_types))
        for (trace_wavelet_idx, trace_type) in zip(flat_patch_wavelet_idxs,
                                                   flat_patch_trace_types)
            if trace_type == AGDLive
                min_wavelet_idx = min(trace_wavelet_idx, min_wavelet_idx)
                max_wavelet_idx = min(trace_wavelet_idx, max_wavelet_idx)
            end
        end
    end
    if (min_wavelet_idx < 0) || (max_wavelet_idx >= length(wavelets))
      error("wavelet_idx must be in [0, length(wavelets))")
    end
end

end
