## Usages

### Generate synthetic structure based on 2D porous medium

step 1: create mask

- command: `createMask2D <binary image data> <nx> <ny>`
- example: `createMask2D origital_structure_small.raw 200 200`

step 2: make synthetic structure

- command: `synthesizeCokedStructure2D <initial patch thickness> <fraction converted to growing material> <growing number>`
- example: `synthesizeCokedStructure2D 2 0.5 5`
