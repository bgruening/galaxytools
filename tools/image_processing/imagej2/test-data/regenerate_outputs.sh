conda_env=mulled-v1-e05ad707e739a59dbca8c6a1fe3f0275ce0ad5649e8bab53f43146c97b48b37e
conda activate $conda_env
# Adjust threshold
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 0.0 129.0 'Default' 'red' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_default.gif' 'gif'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 118.0 255.0 'Percentile' 'over_under' 'no' 'tools/image_processing/imagej2/test-data/blobs_threshold_percentile.gif' 'gif'
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_adjust_threshold_binary_jython_script.py' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 72.0 255.0 'Huang' 'bw' 'yes' 'tools/image_processing/imagej2/test-data/blobs_threshold_huang_dark.gif' 'gif'

# Analyze particles
# Test 1
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Nothing' 'yes' 'yes' 'no' 'no' '' '' 'tools/image_processing/imagej2/test-data/analyze_particles_nothing.tabular'
# Test 2
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Outlines' 'no' 'no' 'no' 'no' 'tools/image_processing/imagej2/test-data/analyze_particles_outlines.gif' 'gif' ''
# Test 3
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' '' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Bare_Outlines' 'no' 'no' 'yes' 'yes' 'tools/image_processing/imagej2/test-data/analyze_particles_bareoutlines.gif' 'gif' ''
# Test 4
ImageJ --ij2 --headless --debug --jython 'tools/image_processing/imagej2/imagej2_analyze_particles_binary_jython_script.py' 'tools/image_processing/imagej2/test-data/analyze_particles_rois.tabular' 'output_log.txt' 'tools/image_processing/imagej2/test-data/blobs.gif' 'no' '0-Infinity' 0.0 1.0 'Overlay_Masks' 'yes' 'yes' 'no' 'no' 'tools/image_processing/imagej2/test-data/analyze_particles_overlaymasks.gif' 'gif' 'tools/image_processing/imagej2/test-data/analyze_particles_nothing.tabular'

