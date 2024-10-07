// colour bed regions according to score - must have columns and score as part of the adapter setup
;(function () {
  class Plugin {
    name = 'BedScorePlugin';
    version = '1.0';

    install(pluginManager) {
      pluginManager.jexl.addFunction('customColor', feature => {
        if (Number(feature.get('score')) > 0) {
          return 'red';
        } else {
          return 'blue';
        }
      })
    }

    configure(pluginManager) {}
  }

  // the plugin will be included in both the main thread and web worker, so
  // install plugin to either window or self (webworker global scope)
  ;(typeof self !== 'undefined' ? self : window).JBrowsePluginBedScorePlugin =
    {
      default: Plugin,
    }
})()
