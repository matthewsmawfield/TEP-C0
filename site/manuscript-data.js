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

function createManuscriptContext() {
    if (typeof require === 'undefined') {
        return {};
    }
    const path = require('path');
    const summaryPath = path.join(__dirname, '..', 'results', 'tep_c0_summary.json');
    const summary = loadJsonIfExists(summaryPath) || {};
    return {
        ...summary.placeholders,
        evidence_gates: summary.evidence_gates || [],
        pipeline: summary.pipeline || {}
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
