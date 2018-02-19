import os	
from pysetta import devel
from pysetta.core.pose import Pose
from pysetta.core.import_pose import pose_from_file

args = "dummy -database "+os.environ['CMAKE_ROSETTA_PATH']+"/database"
devel.init( args.split() )
p = Pose()
p.dump_pdb("test.pdb")
p = pose_from_file("/work/sheffler/1ffw_native.pdb")
p.dump_pdb("test2.pdb")
