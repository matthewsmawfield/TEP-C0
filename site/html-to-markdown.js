#!/usr/bin/env node

const fs = require('fs');
const path = require('path');
const { createManuscriptContext, injectPlaceholders } = require('./manuscript-data.js');

class HTMLToMarkdownConverter {
    decodeEntities(text) {
        return text
            .replace(/&amp;/g, '&')
            .replace(/&lt;/g, '<')
            .replace(/&gt;/g, '>')
            .replace(/&nbsp;/g, ' ')
            .replace(/&ndash;/g, '-')
            .replace(/&mdash;/g, '-')
            .replace(/&Sigma;/g, 'Σ')
            .replace(/&sigma;/g, 'σ')
            .replace(/&oplus;/g, '⊕')
            .replace(/&Epsilon;/g, 'Ε')
            .replace(/&epsilon;/g, 'ε')
            .replace(/&Phi;/g, 'Φ')
            .replace(/&phi;/g, 'φ')
            .replace(/&radic;/g, '√')
            .replace(/&times;/g, '×')
            .replace(/&plusmn;/g, '±')
            .replace(/&infin;/g, '∞')
            .replace(/&alpha;/g, 'α')
            .replace(/&beta;/g, 'β')
            .replace(/&gamma;/g, 'γ')
            .replace(/&delta;/g, 'δ')
            .replace(/&eta;/g, 'η')
            .replace(/&theta;/g, 'θ')
            .replace(/&lambda;/g, 'λ')
            .replace(/&mu;/g, 'μ')
            .replace(/&nu;/g, 'ν')
            .replace(/&pi;/g, 'π')
            .replace(/&rho;/g, 'ρ')
            .replace(/&tau;/g, 'τ')
            .replace(/&chi;/g, 'χ')
            .replace(/&omega;/g, 'ω')
            .replace(/&Omega;/g, 'Ω');
    }

    cleanInline(html) {
        return this.decodeEntities(html)
            .replace(/<(strong|b)[^>]*>([\s\S]*?)<\/(strong|b)>/gi, '**$2**')
            .replace(/<(em|i)[^>]*>([\s\S]*?)<\/(em|i)>/gi, '*$2*')
            .replace(/<code[^>]*>([\s\S]*?)<\/code>/gi, '`$1`')
            .replace(/<sub[^>]*>([\s\S]*?)<\/sub>/gi, '<sub>$1</sub>')
            .replace(/<sup[^>]*>([\s\S]*?)<\/sup>/gi, '<sup>$1</sup>')
            .replace(/<\/?(?!sub|sup)[a-zA-Z][^>]*>/g, '')
            .replace(/\s+/g, ' ')
            .trim();
    }

    tableToMarkdown(tableHtml) {
        const rows = [];
        const rowRegex = /<tr[^>]*>([\s\S]*?)<\/tr>/gi;
        let rowMatch;
        while ((rowMatch = rowRegex.exec(tableHtml))) {
            const cells = [];
            const cellRegex = /<t[hd][^>]*>([\s\S]*?)<\/t[hd]>/gi;
            let cellMatch;
            while ((cellMatch = cellRegex.exec(rowMatch[1]))) {
                cells.push(this.cleanInline(cellMatch[1]).replace(/\|/g, '\\|'));
            }
            if (cells.length) rows.push(cells);
        }
        if (!rows.length) return '';
        const width = rows[0].length;
        const format = (row) => `| ${Array.from({ length: width }, (_, i) => row[i] || '').join(' | ')} |`;
        return `\n\n${format(rows[0])}\n| ${Array.from({ length: width }, () => '---').join(' | ')} |\n${rows.slice(1).map(format).join('\n')}\n\n`;
    }

    htmlToMarkdown(html) {
        return this.decodeEntities(html)
            .replace(/<table[^>]*>[\s\S]*?<\/table>/gi, (match) => this.tableToMarkdown(match))
            .replace(/<h1[^>]*>([\s\S]*?)<\/h1>/gi, '\n# $1\n\n')
            .replace(/<h2[^>]*>([\s\S]*?)<\/h2>/gi, '\n## $1\n\n')
            .replace(/<h3[^>]*>([\s\S]*?)<\/h3>/gi, '\n### $1\n\n')
            .replace(/<h4[^>]*>([\s\S]*?)<\/h4>/gi, '\n#### $1\n\n')
            .replace(/\s*<figcaption[^>]*>([\s\S]*?)<\/figcaption>\s*/gi, '\n\n$1\n\n')
            .replace(/<img[^>]+>/gi, (tag) => {
                const src = (tag.match(/src=["']([^"']+)["']/i) || [])[1];
                const alt = (tag.match(/alt=["']([^"']*)["']/i) || [])[1] || '';
                return src ? `\n\n![${alt}](${src})\n\n` : '';
            })
            .replace(/<blockquote[^>]*>\s*<p[^>]*>([\s\S]*?)<\/p>\s*<\/blockquote>/gi, '\n> $1\n\n')
            .replace(/<p[^>]*>([\s\S]*?)<\/p>/gi, '\n\n$1\n\n')
            .replace(/<li[^>]*>([\s\S]*?)<\/li>/gi, '\n- $1')
            .replace(/<\/?ul[^>]*>/gi, '\n')
            .replace(/<(strong|b)[^>]*>([\s\S]*?)<\/(strong|b)>/gi, '**$2**')
            .replace(/<(em|i)[^>]*>([\s\S]*?)<\/(em|i)>/gi, '*$2*')
            .replace(/<code[^>]*>([\s\S]*?)<\/code>/gi, '`$1`')
            .replace(/<\/?div[^>]*>/gi, '\n')
            .replace(/<\/?section[^>]*>/gi, '\n')
            .replace(/<sub[^>]*>([\s\S]*?)<\/sub>/gi, '<sub>$1</sub>')
            .replace(/<sup[^>]*>([\s\S]*?)<\/sup>/gi, '<sup>$1</sup>')
            .replace(/<\/?(?!sub|sup)[a-zA-Z][^>]*>/g, '')
            .replace(/^[ \t]+/gm, '') // Remove leading indentation
            .replace(/[ \t]+$/gm, '')
            .replace(/\n{3,}/g, '\n\n')
            .trim();
    }

    async convertSiteToMarkdown() {
        const manifestPath = path.join(__dirname, 'manifest.json');
        const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
        const sections = manifest.sections.slice().sort((a, b) => a.order - b.order);
        const context = createManuscriptContext();
        const body = sections
            .map((section) => injectPlaceholders(fs.readFileSync(path.join(__dirname, 'components', section.file), 'utf8'), context))
            .join('\n\n');

        const header = `# ${manifest.title}\n**${manifest.author}**\nVersion: ${manifest.version}\nFirst published: ${manifest.first_published} - Last updated: ${manifest.date}\nDOI: ${manifest.doi}\n\n---\n\n`;
        const markdown = header + this.htmlToMarkdown(body) + '\n';
        const versionOnly = manifest.version.replace(/\s*\([^)]+\)/, '');
        const outputPath = path.join(__dirname, '..', `26-TEP-C0-${versionOnly}-${manifest.codename}.md`);
        fs.writeFileSync(outputPath, markdown, 'utf8');
        console.log(`Markdown saved to: ${outputPath}`);

        // Keep the shared collection archive (../manuscripts/) in sync
        try {
            const sharedArchiveDir = path.join(__dirname, '..', '..', 'manuscripts');
            if (fs.existsSync(sharedArchiveDir)) {
                const archiveName = path.basename(outputPath);
                const paperPrefix = archiveName.split('-')[0];
                for (const staleFile of fs.readdirSync(sharedArchiveDir)) {
                    if (staleFile !== archiveName && staleFile.endsWith('.md') && staleFile.startsWith(`${paperPrefix}-TEP`)) {
                        fs.rmSync(path.join(sharedArchiveDir, staleFile));
                    }
                }
                fs.copyFileSync(outputPath, path.join(sharedArchiveDir, archiveName));
                console.log(`Copied to shared archive: manuscripts/${archiveName}`);
            }
        } catch (archiveError) {
            console.warn(`Could not update shared manuscripts archive: ${archiveError.message}`);
        }
    }
}

if (require.main === module) {
    new HTMLToMarkdownConverter().convertSiteToMarkdown().catch((error) => {
        console.error(error);
        process.exit(1);
    });
}

module.exports = { HTMLToMarkdownConverter };
