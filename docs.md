## Model variants
### Rhodes
The first flavor that we tried in March 2023. LSR neurons and BE IC neurons both send their outputs to separate MOC units. 
Each MOC unit has a lowpass filter and an input-output nonlinearity to generate a gain-control signal on the scale of [0, 1].
The LSR-driven MOC units send those outputs to hair cells over a fairly wide range, whereas the IC-driven MOC units send their signals out to hair cells over a more narrow range. 
At the level of the hair cells, the different MOC inputs are combined multiplicatively to yield a final gain output value in the range [0, 1]. 