# GasNetModel.jl: a Symbolic Problem Builder for Gas Transient Simulations


## Description of the Environment

This package relying on ModelingToolkit.jl allows to compute transient dynamics in a gas network.jl. It has been primarily built to serve as a dependency of GasPwrCoSim.jl co-simulation package but can be used a stand-alone package too.

Data must be loaded to a GasInfo structure, defined as
```
struct GasInfo
    pipes::DataFrame
    nodes::DataFrame
    compressors::DataFrame
    params::DataFrame
end
```
which is then used to populate the model. For more detail on usage, please refer to the tutorial.

## Citation
If you used this package, please  cite our work as
```
@inproceedings{hyett2024differentiable,
  title={Differentiable Simulator For Dynamic \& Stochastic Optimal Gas \& Power Flows},
  author={Hyett, Criston and Pagnier, Laurent and Alisse, Jean and Goldshtein, Igal and Saban, Lilah and Ferrando, Robert and Chertkov, Michael},
  booktitle={2024 IEEE 63rd Conference on Decision and Control (CDC)},
  pages={98--105},
  year={2024},
  organization={IEEE}
}
```






