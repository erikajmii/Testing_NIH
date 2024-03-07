import fs from 'fs';
import { two_nodes } from './two_nodes.mjs';

// Read input parameters from the input file
const inputFilename = 'input.txt';
const outputFilename = 'output.txt';

fs.readFile(inputFilename, 'utf8', async (err, data) => {
    if (err) {
        console.error('Error reading input file:', err);
        return;
    }

    // Split the input data into lines
    const lines = data.trim().split('\n');

    // Process each line
    const results = [];
    for (const line of lines) {
        // Split the line into individual parameters
        const [tdb, rh, met, met2, clo, burn_surface, length_time_simulation] = line.split('\t');

        // Call the function with the provided values
        const t_core = await two_nodes(tdb, rh, met, met2, clo, burn_surface, length_time_simulation);

        // Add the result to the results array
        results.push(t_core);
    }

    // Write the results to the output file
    fs.writeFile(outputFilename, results.join('\n'), 'utf8', (err) => {
        if (err) {
            console.error('Error writing output file:', err);
            return;
        }
        console.log('Results written to', outputFilename);
    });
});
