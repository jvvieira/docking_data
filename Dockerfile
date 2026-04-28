FROM rosettacommons/rosetta:latest

RUN pip install pyrosetta

COPY /src /src

#Create the filtered_data directory if it doesn't exist
RUN mkdir -p /filtered_data

# Create the outputs directory if it doesn't exist
RUN mkdir -p /outputs

USER root

WORKDIR /src