conda_env=mulled-v1-e05ad707e739a59dbca8c6a1fe3f0275ce0ad5649e8bab53f43146c97b48b37e
conda activate $conda_env
# Adjust threshold
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 255.0 'Default' 'red' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_default.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 255.0 'Percentile' 'red' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_percentile.gif' 'gif'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 255.0 'Huang' 'bw' 'yes' 'tools/image_processing/imagej2/test-data/blobs_threshold_huang_dark.gif' 'gif'
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 8.0 255.0 'Manual' 'bw' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_8-255.gif' 'gif'
# Test 5
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 8.0 'Manual' 'bw' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_0-8.gif' 'gif'

# Analyze particles
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Nothing' 'yes' 'yes' 'no' 'no' '' '' 'tools/image_processing/imagej2/test-data/analyze_particles_nothing.tabular'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Outlines' 'no' 'no' 'no' 'no' 'tools/image_processing/imagej2/test-data/analyze_particles_outlines.gif' 'gif' ''
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Bare_Outlines' 'no' 'no' 'yes' 'yes' 'tools/image_processing/imagej2/test-data/analyze_particles_bareoutlines.gif' 'gif' ''
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/analyze_particles_rois.tabular' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Overlay_Masks' 'yes' 'yes' 'no' 'no' 'tools/image_processing/imagej2/test-data/analyze_particles_overlaymasks.gif' 'gif' 'tools/image_processing/imagej2/test-data/analyze_particles_nothing.tabular'

# Skeleton
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/skeletonized_blobs.gif' 'no' 'none' 'no' 'no' 'no' 'tools/image_processing/imagej2/test-data/basic.tabular'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/skeletonized_clown.jpg' 'no' 'shortest_branch' 'no' 'no' 'no' 'tools/image_processing/imagej2/test-data/shortest_branch_basic.tabular'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/skeletonized_blobs.gif' 'no' 'none' 'no' 'yes' 'no' '/tmp/tmpnbw7tzf2/job_working_directory/000/6/outputs/dataset_9ac5bcf4-4e41-4240-871e-7d110c52a7ff.dat'
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_skeleton_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/skeletonized_clown.jpg' 'no' 'shortest_branch' 'no' 'yes' 'no' 'tools/image_processing/imagej2/test-data/shortest_branch_all_yes.tabular'
