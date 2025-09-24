## Build Instructions

Store external dependencies:

```bash
curl -o charts-loader.js https://www.gstatic.com/charts/loader.js
```

Reduce to a single file:

```bash
npm i -g html-inline
html-inline log_template.html -o index.inlined.html
```

Minify file:

```bash
npm i -g html-minifier-terser
html-minifier-terser index.inlined.html \
 --collapse-whitespace \
 --remove-comments \
 --minify-css true \
 --minify-js true \
 -o index.min.html
```
