# Hackathon CVRP

This project was created in the Relaxdays Code Challenge Vol. 1. See https://sites.google.com/relaxdays.de/hackathon-relaxdays/startseite for more information. My participant ID in the challenge was: CC-VOL1-4

## How to run

```bash
git clone https://github.com/Geo5/hackathon-cvrp.git
cd hackathon-cvrp
docker build -t cvrp .
docker run -it cvrp testinstance.json
```

## Algorithmic Idea

The programm solves the CVRP by trasformating it the problem into a linear program and solving it with the free lpsolve software.

The LP formulation used, is based on the second formulation from [this paper](https://hrcak.srce.hr/193543).
