#!/usr/bin/env node

(function () {
    const isNode = typeof module !== 'undefined' && module.exports;

    async function loadBrowserComponents() {
        const sections = document.querySelectorAll('[data-component]');
        for (const section of sections) {
            const file = section.getAttribute('data-component');
            try {
                const response = await fetch(`components/${file}`);
                if (!response.ok) {
                    throw new Error(`${response.status} ${response.statusText}`);
                }
                section.innerHTML = await response.text();
            } catch (error) {
                section.innerHTML = `<p class="component-error">Unable to load ${file}: ${error.message}</p>`;
            }
        }

        if (window.MathJax && window.MathJax.typesetPromise) {
            await window.MathJax.typesetPromise();
        }
    }

    async function buildStaticSite() {
        const fs = require('fs');
        const path = require('path');
        const { HTMLToMarkdownConverter } = require('./html-to-markdown.js');
        const { createManuscriptContext, injectPlaceholders } = require('./manuscript-data.js');

        const siteDir = __dirname;
        const distDir = path.join(siteDir, 'dist');
        const manifest = JSON.parse(fs.readFileSync(path.join(siteDir, 'manifest.json'), 'utf8'));
        let index = fs.readFileSync(path.join(siteDir, 'index.html'), 'utf8');
        const context = createManuscriptContext();

        const components = manifest.sections
            .slice()
            .sort((a, b) => a.order - b.order)
            .map((section) => {
                const componentPath = path.join(siteDir, 'components', section.file);
                const html = fs.existsSync(componentPath)
                    ? injectPlaceholders(fs.readFileSync(componentPath, 'utf8'), context)
                    : `<p>Missing component: ${section.file}</p>`;
                return `<section id="${section.id}" class="manuscript-section" data-section="${section.title}">\n${html}\n</section>`;
            })
            .join('\n\n');

        index = index.replace(
            /<section id="abstract"[\s\S]*?<section id="reproducibility" data-component="8_reproducibility\.html"><\/section>/,
            components
        );
        index = index.replace('<script src="build.js"></script>', '<script>window.__STATIC_MANUSCRIPT__ = true;</script>');

        fs.rmSync(distDir, { recursive: true, force: true });
        fs.mkdirSync(distDir, { recursive: true });
        fs.writeFileSync(path.join(distDir, 'index.html'), index, 'utf8');
        fs.copyFileSync(path.join(siteDir, 'styles.css'), path.join(distDir, 'styles.css'));
        fs.copyFileSync(path.join(siteDir, 'manifest.json'), path.join(distDir, 'manifest.json'));
        
        // Copy public directory to dist
        const publicDir = path.join(siteDir, 'public');
        const publicDest = path.join(distDir, 'public');
        if (fs.existsSync(publicDir)) {
            copyRecursiveSync(publicDir, publicDest, fs, path);
        }
        
        copyResults(path.join(siteDir, '..', 'results'), distDir, fs, path);
        copyDataArtifacts(path.join(siteDir, '..'), distDir, fs, path);

        for (const file of ['CITATION.cff', 'CITATION.bib', 'README.md', 'LICENSE', 'VERSION.json']) {
            const src = path.join(siteDir, '..', file);
            if (fs.existsSync(src)) {
                fs.copyFileSync(src, path.join(distDir, file));
            }
        }

        const converter = new HTMLToMarkdownConverter();
        await converter.convertSiteToMarkdown();

        console.log(`Built TEP-C0 static manuscript: ${path.join(distDir, 'index.html')}`);
    }

    function copyResults(resultsDir, distDir, fs, path) {
        if (!fs.existsSync(resultsDir)) {
            return;
        }
        const figuresSrc = path.join(resultsDir, 'figures');
        const figuresDest = path.join(distDir, 'figures');
        if (fs.existsSync(figuresSrc)) {
            fs.mkdirSync(figuresDest, { recursive: true });
            for (const file of fs.readdirSync(figuresSrc)) {
                if (file.endsWith('.png')) {
                    fs.copyFileSync(path.join(figuresSrc, file), path.join(figuresDest, file));
                }
            }
        }

        const resultsDest = path.join(distDir, 'results');
        fs.mkdirSync(resultsDest, { recursive: true });
        for (const file of fs.readdirSync(resultsDir)) {
            if (file.endsWith('.json') || file.endsWith('.csv')) {
                fs.copyFileSync(path.join(resultsDir, file), path.join(resultsDest, file));
            }
        }
    }

    function copyDataArtifacts(projectDir, distDir, fs, path) {
        const resultsDir = path.join(projectDir, 'results');
        const payloadFiles = [
            path.join(resultsDir, 'step_000_data_download.json'),
            path.join(resultsDir, 'step_000_data_ingestion.json'),
        ];
        const artifactPaths = new Set();

        for (const payloadFile of payloadFiles) {
            if (!fs.existsSync(payloadFile)) {
                continue;
            }
            const payload = JSON.parse(fs.readFileSync(payloadFile, 'utf8'));
            for (const artifactPath of Object.values(payload.artifacts || {})) {
                if (typeof artifactPath === 'string' && artifactPath.startsWith('data/')) {
                    artifactPaths.add(artifactPath);
                }
            }
        }

        for (const artifactPath of artifactPaths) {
            const src = path.join(projectDir, artifactPath);
            if (!fs.existsSync(src) || !fs.statSync(src).isFile()) {
                continue;
            }
            const dest = path.join(distDir, artifactPath);
            fs.mkdirSync(path.dirname(dest), { recursive: true });
            fs.copyFileSync(src, dest);
        }
    }

    function copyRecursiveSync(src, dest, fs, path) {
        const stat = fs.statSync(src);
        if (stat.isDirectory()) {
            if (!fs.existsSync(dest)) {
                fs.mkdirSync(dest, { recursive: true });
            }
            const entries = fs.readdirSync(src);
            for (const entry of entries) {
                copyRecursiveSync(path.join(src, entry), path.join(dest, entry), fs, path);
            }
        } else {
            fs.copyFileSync(src, dest);
        }
    }

    if (isNode) {
        module.exports = { buildStaticSite };
        if (require.main === module) {
            buildStaticSite().catch((error) => {
                console.error(error);
                process.exit(1);
            });
        }
    } else if (!window.__STATIC_MANUSCRIPT__) {
        document.addEventListener('DOMContentLoaded', loadBrowserComponents);
    }
})();
