For this exercise we will be running a truncated version of ClimaAtmos.jl, the atmospheric component model of the CliMA ESM introduced in the last lecture

# Setup
- Before/at the beginning of the lecture, it would be helpful if you could
    - Download the latest Julia Version: https://julialang.org/downloads/
        - tip: make Julia callable from the Terminal ln -s ../$APPLICATION_PATH/julia /usr/local/bin/julia
    - Create a GitHub account (if you don’t already have one)
    - For code modification and figure display we’ll be using VSCode, but you can use your favorite software

# Initial exploration
- Git clone https://github.com/LenkaNovak/ClimaAtmos_ESE136.jl
    - Explore the structure of the model
- Run the model for a short time (e.g., 2 days) to ensure it’s setup correctly
    - `cd examples`
    - `julia --project —threads 8`
    - `]instantiate`
    - `include("hybrid/driver.jl")`
    - Explore output using the post processing scripts called at the end of driver.jl (or you can only apply the regridding, and use your own code to read and explore the NetCDF files)
- let's make a longer run longer (e.g., 100 days - this may take ~1h with 1 CPU process with 8 threads), to get behavior beyond the spin-up period
    - In the driver
    - Set t_end to 100days
    - Set dt_save_to_disk to 5days
    - If using the plot_time_mean for post processing, change the day_range to collect(50:1:100)
    - Inspect the time-average (”climatology”) plots. How does this compare to Earth’s climatology?

# Experiment 1: Effect of rotation
- Rotation change: Change the rotation of the Earth (`Omega` in examples/hybrid/parameter_set.jl) both increasing and decreasing it (e.g. by a factor of 2). Does this agree with what you would expect, based on existing literature? How does this compare to the energy balance model results from the last lab?

# Experiment 2: Hypo-hydrostatic scaling
- Many aspects of the large-scale circulation can be reproduced on a smaller planet (if the rotation and frictional timescales are adapted accordingly). This allows cheaper simulations, and/or better resolved simulations which can coarsely resolve convection. Can you reproduce the Earth-sized aquaplanet climate on a dwarf planet?
    - Lets use a scaling factor of `X`= 8
    - radius / `X` (`planet_radius` in examples/hybrid/parameter_set.jl)
    - corollas parameter x `X` (`Omega` in examples/hybrid/parameter_set.jl)
    - all inverse damping timescales x `X` (`parsed_args["C_E"]` in examples/hybrid/parameter_set.jl)
    - Since we want to damp smaller wavelengths by hyper diffusion, we reduce its coefficient  (`parsed_args["kappa_4"]`  in driver.jl ) to 2e14 m4/s
    - You will also need to reduce the time step (try 40s) and make the simulation shorter (`t_end` = 10days)
- based on your understanding of the literature, explain why this setup should yield a similar climate to the large planet, and when this scaling breaks down.

# Experiment 3 (optional): Changing the surface temperature
- Changing the `T_sfc` in the default_cache() - create a localized Gaussian heating patch in the tropics. What can you see?

# Troubleshooting
- if the above set up doesn't work (e.g., we've had some issues with WindowsOS), there are three other ways:
    - cloning and running this tutorial on the HPC cluster (ask if you don't have access yet)
        - careful not to use the login node. Please request a compute node before calculating anything (e.g., `srun --pty -t 2:00:00 -n 1 /bin/bash -l` for 2 hours on one process)
    - cloning this tutorial and running the examples using a SLURM script - we recommend parallelizing using multithreading, but you can try out parallelizing using MPI too (see `examples/README.md` for instructions).
    - using [this GoogleColab](https://colab.research.google.com/drive/1d6APTLOeAEoal67nT1ZqGErZRUyCeb1S?usp=sharing) (though this should be your last resort as the full instantiation has to be done with each new run time).
        - remember to copy this to your drive, and refresh after running the Julia install cell
        - to change parameters of the cloned `ClimaAtmos_ESE136.jl` (see below), you'll need to access this in your drive (by default it should be saved in your `Colab Notebooks` directory).

# References
- Aquaplanet experiments
    - https://www.gfdl.noaa.gov/bibliography/related_files/dmwf0601.pdf
- Rotation change
    - https://www.gfdl.noaa.gov/bibliography/related_files/ann0201.pdf
    - https://journals.ametsoc.org/view/journals/clim/34/9/JCLI-D-20-0533.1.xml
- Hypo-hydrostatic scaling
    - https://www.atmos.washington.edu/~dargan/papers/gfhpv07.pdf
    - https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022MS003527