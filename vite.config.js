// vite.config.js
import { defineConfig } from "vite";
import anywidget from "@anywidget/vite";
import license from 'rollup-plugin-license';

export default defineConfig({
	build: {
		outDir: "./src/cobramod/static",
		lib: {
			entry: {
                escher: "./js/escher.ts",
                force_graph: "./js/force_graph.ts"
            },
			formats: ["es"],
		},
	},
	plugins: [
		anywidget(),
		license({
		    sourcemap: true,
			thirdParty: {
			  output: './src/cobramod/static/LICENSES.txt',
			  allow: '(MIT OR ISC OR BSD-3-Clause)',
			},
	    }),

	],
});