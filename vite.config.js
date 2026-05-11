// vite.config.js
import { defineConfig } from "vite";
import anywidget from "@anywidget/vite";
import license from 'rollup-plugin-license';
import { build as esbuildBundle } from 'esbuild';
import { resolve } from 'path';

// flux_network must be built as a single self-contained file — anywidget loads it as
// a Blob URL, so relative chunk imports produced by Rollup/Vite would fail.
// Build it with esbuild (which never splits chunks) after the main Vite build.
function bundleFluxNetwork() {
	return {
		name: 'bundle-flux-network',
		apply: 'build',
		async closeBundle() {
			await esbuildBundle({
				entryPoints: [resolve('./js/flux_network/index.ts')],
				bundle: true,
				format: 'esm',
				outfile: resolve('./src/cobramod/static/flux_network.mjs'),
			});
		}
	};
}

export default defineConfig({
	resolve: {
		dedupe: ['three'],
	},
	build: {
		outDir: "./src/cobramod/static",
		lib: {
			entry: {
                escher: "./js/escher.ts",
                force_graph: "./js/force_graph.ts",
                // flux_network is built by the bundleFluxNetwork plugin (esbuild)
            },
			formats: ["es"],
		},
	},
	plugins: [
		anywidget(),
		bundleFluxNetwork(),
		license({
		    sourcemap: true,
			thirdParty: {
			  output: './src/cobramod/static/LICENSES.txt',
			  allow: '(MIT OR ISC OR BSD-3-Clause)',
			},
	    }),
	],
});
