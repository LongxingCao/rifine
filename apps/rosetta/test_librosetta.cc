#include <iostream>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/annotated_sequence.hh>

int main(int argc, char *argv[])
{
	::devel::init( argc, argv );
	std::cout << "hello rosetta" << std::endl;

	core::pose::Pose pose;
	core::pose::make_pose_from_sequence(pose,"THANKSEVAN",core::chemical::FA_STANDARD);
	pose.dump_pdb("test.pdb");

	return 0;
}