import { two_nodes } from './two_nodes.mjs';

// Call the function with the provided values
//environemal temp, rh, total heat productin, relative heat production, burn surface, duration
const t_core = await two_nodes(39.079, 21.20, 530, 4.90, 0.15, 40.00, 60);

// Log the result to the console
console.log("t_core:", t_core);
