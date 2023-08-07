The **Statistics View** provides more details about the ligands of one selected cluster or the whole docking, depending on if a selection in the **Docking Overview** is made or not. Here the cluster-specific information is presented as a table enabling to explore the **_ligand statistics_**, **_scores_**.

The table itself is scrollable and includes values like _min_ and _max score_ that were found for all ligand poses which are included in the selected cluster and counted to one specific ligand. Statistical values are predominantly averaged values for the selected cluster or averaged over all clusters where poses of that ligand are present.

- ### Columns & Sorting

  The heading of the table will always show which cluster is currently selected and an **on-hover tooltip** is provided for each of the column names with a detailed description. Any column with alphanumerical **values can be sorted**. It offers the identification of most scoring ligands. While hovering over the column names an arrow appears indicating that it is sortable. A click activates the sorting and can also invert the ordering (ascending/descending).

- ### Special Columns

  The rows are selectable by a click and a selected ligand is highlighted by the change of the row’s background color (🟩). Besides columns showing values, there are further, e.g. the **second **column**, containing the names of the ligands, which work simultaneously as **web-link to the ZINC database** pages of the individual ligands. From this web page, more detailed information can be obtained. As a complement to this and to help the user to get an impression of **molecular structure** the **3rd column** provides a two-dimensional representation of the ligands. The structure is displayed as a structural formula using the Natta projection to show the spatial orientation of chemical groups. Because some ligands can be larger and hard to read, a **magnifying tooltip\*\* is implemented. When the user points with their cursor to a structural formula a magnified version appears on the left side of it.

- ### Scores

  Moreover, the table contains a colored, rounded background for the mentioned _max_ and _min scores_. It indicates directly if the **AutoDock-score/free binding energy** is high (more negative) or low (less negative). We applied the _inferno color map_ to it going from min values (almost black) to max values (very light yellow). While hovering a score a color legend is shown with the whole color map and corresponding score range together with a green indicator arrow pointing to the exact value. This legend also shows up by mouseover of the blue **LEGEND** button at the top right of the table.

- ### Search & Filter

  Besides this basic option to rearrange the table content by sorting after column values, we extended the **Statistics View** with **search** and **filter** possibilities. The first field above the table **Filter by Interaction Type** allows the user to filter by a specific interaction type. Further, the user can refine the filtering or just filter the table by the second field **Search Table** to the left of the previously described field. The search field takes any alphanumeric input, which makes it possible, e.g., to search for a specific ZINC name or certain scores.
