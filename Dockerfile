FROM continuumio/miniconda3

WORKDIR /cvrp

# Create the environment:
COPY . .
RUN conda create -n cvrp -y

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "cvrp", "/bin/bash", "-c"]

# Make sure the environment is activated:
RUN conda install -c conda-forge --file requirements.txt -y

# The code to run when container is started:
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "cvrp", "python", "cvrp.py"]
