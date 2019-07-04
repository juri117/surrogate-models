# Surrogate Models

In this repository I publish the python code, that was part of my master thesis. The thesis can be found [here](docs/Masterarbeit.pdf), however its in German though, sry. :/

## content

### implementations of surrogate models

The folder [mylibs](mylibs) contains implementations of surrogate models. The following methods are used:

* polynomial
* radial basis function
* kriging

moreover some methods for generating sample points are implemented:

* structured sample
* halton
* latin hyper cube

### Examples

The folder [examples](examples) shows how the above libraries can be used both in 2D and 3D spaces.

### wing-construction

The construction of a wing was the sample task I constructed for testing and comparing the surrogate methods. A wing-box was meant to be optimized for its weight, by varying its numbers of ribs and shell thickness. To make sure the structure is strong enough a FEM analysis is used (both Calculix and Abaqus APIs are implemented to do this).