from shutil import rmtree
import new_Constants

import os
import sys


def get_clusters_and_centers(seq_file_path, is_precluster_mode = False):
    cluster_ids_to_centers_and_cluster_seqs = dict()
    # file locations
    main_dir_path = os.path.dirname(os.path.realpath(__file__))


    # SET UP CONFIGURATIONS

    try:
        print("hi")

    except KeyboardInterrupt:
        print()
        print('Process aborted due to keyboard interrupt')
    except SystemExit as sys_exit:
        if sys_exit.code != 0:
            print(sys_exit.code)
    except:
        print()
        print('Process aborted due to error occurred: {}'.format(sys.exc_info()[1]))
    # finally:
    #     Precluster.clear_temp_data()