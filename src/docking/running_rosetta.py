import pyrosetta
import os

pyrosetta.init('-mute all')  # Mute Rosetta output for cleaner logs

def main():
    print("Running Rosetta...")
    files = os.listdir("/filtered_data")
    print(f"Files in /filtered_data: {len(files)}")
    
    # Load Pose
    pose = pyrosetta.pose_from_pdb(f"/filtered_data/{files[0]}") 

    print(pose.get_chains())  # Print variant types to check if the pose is loaded correctly
    # print(pose.split_by_chain(1).sequence())  # Example: split chain 1
    # print(pose.split_by_chain(2).sequence())  # Example: split chain 2

if __name__ == "__main__":
    main()