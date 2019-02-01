require(['mathjax'], function(MathJax) {
  MathJax.Hub.Config({
    TeX: {
      Macros: {
        RR: "{\\bf R}",
        bold: ["{\\bf #1}",1]
      }
    }
  });

  //MathJax.Ajax.loadComplete("[MathJax]/config/local/mathjax.js");
})
