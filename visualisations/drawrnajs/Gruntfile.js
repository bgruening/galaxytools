// NOTE: use 'sudo npm install .', then 'grunt' to use this file
module.exports = function(grunt)
{
    grunt.initConfig({
        pkg: grunt.file.readJSON('package.json'),
        
          uglify: {
            options: {
                    mangle: true,
                    compress: true
            },
            dist: {
                files: {
                    'static/function.os_to_dotbracket.min.js': 'src/function.os_to_dotbracket.js',
                    'static/function.parse_ps.min.js': 'src/function.parse_ps.js',
                    'static/function.redraw.min.js': 'src/function.redraw.js',
                    'static/function.draw_sequence_selector.min.js': 'src/function.draw_sequence_selector.js',
                    'static/function.get_fasta_index.min.js': 'src/function.get_fasta_index.js',
                    'static/function.get_dotbracket_index.min.js': 'src/function.get_dotbracket_index.js',
                    'static/function.get_connectivitytable_index.min.js': 'src/function.get_connectivitytable_index.js',
                    'static/drawrnajs/drawrnajs.0.3.5.min.js': 'src/drawrnajs/drawrnajs.0.3.5.js'
                }
            }
        },
    });
    
    // Load the plugin that provides the "uglify" task.
    grunt.loadNpmTasks('grunt-contrib-uglify');
    
    // Default task(s).
    grunt.registerTask('default', ['uglify']);
};
