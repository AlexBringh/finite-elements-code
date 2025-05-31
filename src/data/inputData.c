#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "inputData.h"

int readCSV(const char *filename, node **nodes_out, int *nnodes_out, int *nnodesElement, int *nelements, int *gp, material **materials_out, int *nmaterials_out, int **uFixed_out, double **fApplied_out) 
{

    FILE *fp = fopen(filename, "r");
    if (!fp) 
    {
        perror("Failed to open CSV file");
        return -1;
    }

    char line[256];
    int section = 0;
    int count_nodes = 0, count_mats = 0;
    int nnodes = 0;
    int dof = 0;
    int system_values[6] = {0}; // {nelements, nnodes, nnodeselement, gp, dof, dim}

    node *nodes = NULL;
    material *materials = malloc(100 * sizeof(material));
    if (!materials) 
    {
        perror("malloc failed for materials");
        fclose(fp);
        return -1;
    }

    while (fgets(line, sizeof(line), fp)) 
    {
        line[strcspn(line, "\n")] = 0;
        if (strlen(line) == 0) continue;

        if (strcmp(line, "-1") == 0) 
        {
            section++;
            continue;
        }

        // Section 0: Material properties
        if (section == 0) 
        {
            material m;
            int scanned = sscanf(line, "%d,%lf,%lf,%lf,%lf", &m.id, &m.E, &m.v, &m.H, &m.yield);
            if (scanned == 5) 
            {
                materials[count_mats++] = m;
            }
        }

        // Section 1: System info
        else if (section == 1) 
        {
            sscanf(line, "%d,%d,%d,%d,%d,%d",
                   &system_values[0], &system_values[1], &system_values[2],
                   &system_values[3], &system_values[4], &system_values[5]);
            nnodes = system_values[1];
            dof    = system_values[4];
            *gp = system_values[3];
            *nelements = system_values[0];
            *nnodesElement = system_values[2];

            // Allocate displacement and load vectors
            *uFixed_out   = calloc(nnodes * dof, sizeof(double));
            *fApplied_out = calloc(nnodes * dof, sizeof(double));
            if (!*uFixed_out || !*fApplied_out) 
            {
                perror("calloc failed for uFixed or fApplied");
                fclose(fp);
                free(materials);
                return -1;
            }
        }

        // Section 2: Node data
        else if (section == 2) 
        {
            if (nodes == NULL) 
            {
                nodes = malloc(nnodes * sizeof(node));
                if (!nodes) 
                {
                    perror("malloc failed for nodes");
                    fclose(fp);
                    free(materials);
                    return -1;
                }
            }

            if (count_nodes >= nnodes) 
            {
                fprintf(stderr, "Warning: more nodes than declared (%d)\n", nnodes);
                break;
            }

            node n;
            int scanned = sscanf(line, "%d,%lf,%lf,%d", &n.id, &n.x, &n.y, &n.material);
            if (scanned == 4) 
            {
                n.dof = system_values[4];
                n.dim = system_values[5];
                nodes[count_nodes++] = n;
            }
        }

        // Section 4: Fixed Displacements
        else if (section == 4) 
        {
            int nodeid;
            int ux, uy;
            int scanned = sscanf(line, "%d,%d,%d", &nodeid, &ux, &uy);
            if (scanned == 3) 
            {
                (*uFixed_out)[nodeid * dof + 0] = ux;
                (*uFixed_out)[nodeid * dof + 1] = uy;
            }
        }

        // Section 5: Applied Loads
        else if (section == 5) 
        {
            int nodeid;
            double fx, fy;
            int scanned = sscanf(line, "%d,%lf,%lf", &nodeid, &fx, &fy);
            if (scanned == 3) 
            {
                (*fApplied_out)[nodeid * dof + 0] = fx;
                (*fApplied_out)[nodeid * dof + 1] = fy;
            }
        }

        // Stop if we're past Section 5
        else if (section > 5) 
        {
            break;
        }
    }

    fclose(fp);
    *nodes_out = nodes;
    *nnodes_out = count_nodes;
    *materials_out = materials;
    *nmaterials_out = count_mats;

    return 0;
}


int readElements(const char *filename, quadElement **elements_out, int *nelements_out, node *nodes, int nnodes, int gp, int nnodes_per_elem, int dof) 
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Failed to open file in readElements");
        return -1;
    }

    char line[256];
    int section = 0;
    int element_count = 0;
    int capacity = 100;

    quadElement *elements = malloc(capacity * sizeof(quadElement));
    if (!elements) 
    {
        perror("malloc failed for elements");
        return -1;
    }

    // Reset file pointer and re-parse to section 3
    rewind(fp);
    while (fgets(line, sizeof(line), fp)) 
    {
        line[strcspn(line, "\n")] = 0;
        if (strcmp(line, "-1") == 0) 
        {
            section++;
            continue;
        }

        // Element connectivity section
        if (section == 3) 
        { 
            if (element_count >= capacity) 
            {
                capacity *= 2;
                elements = realloc(elements, capacity * sizeof(quadElement));
                if (!elements) {
                    perror("realloc failed for elements");
                    return -1;
                }
            }

            quadElement e;
            e.gp = gp;
            e.nnodes = nnodes_per_elem;
            e.dof = dof;

            // Allocate and parse node IDs
            e.nodeids = malloc(nnodes_per_elem * sizeof(int));
            if (!e.nodeids) {
                perror("malloc failed for nodeids");
                return -1;
            }

            int id, n0, n1, n2, n3;
            int scanned = sscanf(line, "%d,%d,%d,%d,%d", &id, &n0, &n1, &n2, &n3);
            if (scanned != 5) continue;

            e.id = id;
            e.nodeids[0] = n0;
            e.nodeids[1] = n1;
            e.nodeids[2] = n2;
            e.nodeids[3] = n3;

            // Allocate and fill coords (flattened)
            e.coords = malloc(nnodes_per_elem * dof * sizeof(double));
            if (!e.coords) {
                perror("malloc failed for coords");
                return -1;
            }

            for (int i = 0; i < nnodes_per_elem; i++) {
                int node_id = e.nodeids[i];
                int found = 0;
                for (int j = 0; j < nnodes; j++) {
                    if (nodes[j].id == node_id) {
                        e.coords[i * 2 + 0] = nodes[j].x;
                        e.coords[i * 2 + 1] = nodes[j].y;
                        found = 1;
                        break;
                    }
                }
                if (!found) {
                    fprintf(stderr, "Node ID %d not found for element %d\n", node_id, id);
                    return -1;
                }
            }

            int stressStrainSize = gp * 3;

            // Allocate internal stress/strain arrays
            e.sigma             = calloc(stressStrainSize, sizeof(double));
            e.epsilonP          = calloc(stressStrainSize, sizeof(double));
            e.epsilonBarP       = calloc(gp, sizeof(double));
            e.trialSigma        = calloc(stressStrainSize, sizeof(double));
            e.trialEpsilonP     = calloc(stressStrainSize, sizeof(double));
            e.trialEpsilonBarP  = calloc(gp, sizeof(double));

            elements[element_count] = e;
            element_count += 1;
        }

        if (section > 3) break; // done reading element data
    }

    *elements_out = elements;
    *nelements_out = element_count;
    return 0;
}

