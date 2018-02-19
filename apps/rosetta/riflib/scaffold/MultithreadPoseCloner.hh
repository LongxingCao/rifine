// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#ifndef INCLUDED_riflib_scaffold_MultithreadPoseCloner_hh
#define INCLUDED_riflib_scaffold_MultithreadPoseCloner_hh


#include <scheme/types.hh>

#include <core/pose/Pose.hh>
#include <mutex>
#include <thread>



namespace devel {
namespace scheme {

// This class allows one to make the minimum copies of a pose while
//  in a multithreaded environment where all threads want to clone the
//  pose. 
//
// Works as follows:
//  If one of the poses is available. Lock and clone it
//  If none of the poses are available. Wait for a pose, duplicate it, 
//     store that pose in poses. And then retry again.

// Theoretically the max number of poses == number of threads
// Alternatively, it's possible the list will never grow past 1

struct MultithreadPoseCloner {

    std::vector<core::pose::PoseCOP> poses_;
    std::deque<std::mutex> pose_mutexes_;

    std::mutex vector_mutex_;

    MultithreadPoseCloner() : need_more_poses_(0) {}

    MultithreadPoseCloner(core::pose::PoseCOP pose) : need_more_poses_(0) {
        add_pose( pose );
    }

    void
    add_pose( core::pose::PoseCOP pose ) {
        std::lock_guard<std::mutex> guard( vector_mutex_ );
        poses_.push_back(pose);
        pose_mutexes_.resize( poses_.size() );
    }

    core::pose::PoseCOP
    get_pose() {
        runtime_assert( poses_.size() > 0);

        core::pose::PoseCOP to_return;

        while ( ! to_return ) {
            int got_lock_at = -1;
            {
                std::lock_guard<std::mutex> guard( vector_mutex_ );
                for ( int i = 0; i < pose_mutexes_.size(); i++ ) {
                    if ( pose_mutexes_[i].try_lock() ) {
                        got_lock_at = i;
                        break;
                    }
                }
            }
            if (got_lock_at != -1) {
                to_return = clone_a_pose( poses_[got_lock_at] );
                pose_mutexes_[got_lock_at].unlock();

                bool should_duplicate = false;
                {
                    std::lock_guard<std::mutex> guard( lock_need_more_poses_ );
                    if (need_more_poses_ > 0) {
                        need_more_poses_ -= 1;
                        should_duplicate = true;
                    }  
                }
                if (should_duplicate) {
                    core::pose::PoseCOP new_pose = clone_a_pose( to_return );
                }
            } else {
                std::this_thread::sleep_for( std::chrono::milliseconds {1});
            }
        }

        return to_return;

    }

    // Don't be super smart/clever here. Just wait for a mutex on the last
    // element
    void
    duplicate_a_pose() {
        int pose_to_lock = 0;
        {
            std::lock_guard<std::mutex> guard( vector_mutex_ );
            pose_to_lock = pose_mutexes_.size() - 1;
        }
        {
            std::lock_guard<std::mutex> guard( lock_need_more_poses_ );
            need_more_poses_ += 1;
        }
        pose_mutexes_[pose_to_lock].lock();
        core::pose::PoseCOP new_pose = clone_a_pose( poses_[pose_to_lock] );
        pose_mutexes_[pose_to_lock].unlock();

        add_pose( new_pose );
    }

    uint64_t
    size() {
        std::lock_guard<std::mutex> guard( vector_mutex_ );
        return poses_.size();
    }

    static
    core::pose::PoseCOP
    clone_a_pose( core::pose::PoseCOP pose ) {
        core::pose::PoseOP to_return = make_shared<core::pose::Pose>();
        to_return->detached_copy( *pose );
        return to_return;
    }


private:

    std::mutex lock_need_more_poses_;
    int need_more_poses_;
};


}}



#endif