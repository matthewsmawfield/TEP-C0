function loadJsonIfExists(path) {
    if (typeof require === 'undefined') {
        return null;
    }
    const fs = require('fs');
    if (!fs.existsSync(path)) {
        return null;
    }
    return JSON.parse(fs.readFileSync(path, 'utf8'));
}

function loadSharedSnippet(path) {
    if (typeof require === 'undefined') {
        return null;
    }
    const fs = require('fs');
    if (!fs.existsSync(path)) {
        return null;
    }
    return fs.readFileSync(path, 'utf8').trim();
}

function createManuscriptContext() {
    if (typeof require === 'undefined') {
        return {};
    }
    const path = require('path');
    const summaryPath = path.join(__dirname, '..', 'results', 'tep_c0_summary.json');
    const summary = loadJsonIfExists(summaryPath) || {};

    // Load shared corpus snippets from core/
    const screeningNoticePath = path.join(__dirname, '..', 'core', 'screening_projection_notice.html');
    const screeningNotice = loadSharedSnippet(screeningNoticePath);

    return {
        ...summary.placeholders,
        evidence_gates: summary.evidence_gates || [],
        pipeline: summary.pipeline || {},
        SCREENING_PROJECTION_NOTICE: screeningNotice || '<!-- screening_projection_notice.html not found -->'
    };
}

function injectPlaceholders(template, context = createManuscriptContext()) {
    return template.replace(/\{\{\s*([A-Za-z0-9_]+)\s*\}\}/g, (match, key) => {
        const value = context[key];
        if (value === undefined || value === null) {
            return match;
        }
        return String(value);
    });
}

if (typeof module !== 'undefined') {
    module.exports = { createManuscriptContext, injectPlaceholders };
}
