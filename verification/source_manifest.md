# Source manifest

## Frozen analytical source

| Item | Retrieval/freeze date | Public location | SHA-256 |
| --- | --- | --- | --- |
| OWID Energy dataset (`owid-energy-data.csv`) | 1 September 2026 | <https://owid-public.owid.io/data/energy/owid-energy-data.csv> | `77b3db513f02f5fffb69fe02832907ce70b01d3906fc2c5dd40fa47e3ee7d0f3` |
| OWID Energy codebook (`owid-energy-codebook.csv`) | 1 September 2026 | <https://github.com/owid/energy-data/blob/master/owid-energy-codebook.csv> | `3cc9b7db0d921496e2988568ce3aee5ed41f50431dd234a0663b5f0a4b2e32bb` |

Both hashes were independently reverified on 2 September 2026 before this
reproducibility archive was assembled.

## Provider lineage used in the study

The frozen OWID codebook identifies the electricity series used in the study as
harmonized from Ember's Yearly Electricity Data and the Energy Institute
Statistical Review of World Energy. The codebook remains the authoritative
field-level record of definitions, units, processing notes, and sources.

## Analytical variables

| Analytical role | OWID variable or construction | Unit |
| --- | --- | --- |
| Renewable generation | `renewables_electricity` | TWh |
| Fossil generation | `fossil_electricity` | TWh |
| Total generation | `electricity_generation` | TWh |
| Electricity demand | `electricity_demand` | TWh |
| Fossil generation share | `fossil_share_elec` | % |
| Net-import share | `net_elec_imports_share_demand` | % of demand |
| Wind + solar share | (`wind_electricity` + `solar_electricity`) / total generation | % |
| Nuclear share | `nuclear_electricity` / total generation | % |

## Population and period

The population comprises the 19 country members of the G20. The analytical
period is 2000–2024; 1999 is used only as a lag buffer. The European Union and
African Union are supranational members and are not country-level observational
units. The balanced analytical panel contains 475 country-years.
