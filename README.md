# VirusGameRust
[VirusGame](https://github.com/carykh/VirusGame/) with a spin on it, rewritten in Rust and Bevy.

<img src="https://user-images.githubusercontent.com/11985242/175776480-ff635e2d-1c76-466f-9226-ea89f5a03cb3.png" width="612" height="540">

## Todos

- [x] - Implement cells/walls that do nothing lol
- [x] - Implement particle functionality
- [x] - Implement food spawner
- [x] - Implement simple food -> energy -> waste cycle with a spin on it (eat -> get energy -> consume energy -> free itself up from waste to get more energy)
- [x] - Implement simple genes (OG VirusGame basis) functionality
- [ ] - Implement cell walls [**BLOCKED**](https://github.com/bevyengine/bevy/issues/5081)
- [ ] - Mutations [*Partially*]
- [X] - Implement codons visual representation
- [ ] - Implement UGOs
- [ ] - Implement more advanced waste/ugo removal system
- [ ] - Implement advanced genes functionality (Change gene reading direction)
- [ ] - Implement multiplication for cells
- [ ] - GUI editor for cell genes/ugos
- [ ] - Speed controls (play/pause, faster, slower)
- [ ] - Simulation settings in general
- [ ] - Stats (alive, mutations, "natural" vs "replicated")

## Nice to haves

- [ ] - Write up what genes do
- [ ] - Graphs for stats
- [ ] - Gene stats + graphs (most common snips, gene strings)
- [ ] - More complex genes (regulate mutation rate, cell ticking rate etc.)
- [ ] - Food production/storage by cells themselves
- [ ] - Toggle to use only original genes + functionality (maybe)
- [ ] - World exporting/importing
- [ ] - Cell/organism importing/exporting
- [ ] - Cell communication (hormones/pheromones)



## Run instructions
1. Get Rust
2. For quicker builds get [LLD](https://github.com/carykh/VirusGame/)
3. `cargo run` for developing, `cargo run --release` to run with all optimizations.