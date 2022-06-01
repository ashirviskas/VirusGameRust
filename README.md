# VirusGameRust
[VirusGame](https://github.com/carykh/VirusGame/) with a spin on it, rewritten in Rust and Bevy.

## Todos
[ ] - Implement cells/walls that do nothing lol
[ ] - Implement particle functionality
[ ] - Implement simple food -> energy -> waste cycle with a spin on it (eat -> get enerty -> consume energy -> to get more energy, free itself up from waste)
[ ] - Implement simple genes (OG VirusGame basis) functionality
[ ] - Implement cell walls
[ ] - Mutations
[ ] - Implement codons visual representation
[ ] - Implement UGOs
[ ] - Implement more advanced waste/ugo removal system
[ ] - Implement advanced genes functionality ()
[ ] - Implement multiplication for cells
[ ] - GUI editor for cell genes/ugos
[ ] - Speed controls (play/pause, faster, slower)
[ ] - Simulation settings in general
[ ] - Stats (alive, mutations, "natural" vs "replicated")
## Nice to haves
[ ] - Write up what genes do
[ ] - Graphs for stats
[ ] - Gene stats + graphs (most common snips, gene strings)
[ ] - More complex genes (regulate mutation rate, cell ticking rate, etc.)
[ ] - Food production/storage by cells themselves
[ ] - Toggle to use only original genes + functionality (maybe)
[ ] - World exporting/importing
[ ] - Cell/organism importing/exporting
[ ] - Cell communication (hormones/pheromones)



## Run instructions
1. Get Rust
2. For quicker builds get [LLD](https://github.com/carykh/VirusGame/)
3. `cargo run` for developing, `cargo run --release` to run with all optimizations.