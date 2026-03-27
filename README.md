# spatial-transcriptomics-scripts
Ongkeko Lab spatial transcriptomics data processing pipeline

in st_master, define path and name variables

flow chart of script dependencies
<img width="1126" height="768" alt="image" src="https://github.com/user-attachments/assets/b6bd4170-672e-414d-8820-f13739305e72" />
https://lucid.app/lucidchart/9ca23dbf-cc43-4b66-b452-ce79e9b11e3e/edit?viewport_loc=-547%2C-325%2C3078%2C1476%2C0_0&invitationId=inv_80ce0e20-4ee8-4184-8f4a-6fc41391dc04

Make sure each sample's file structure is organized as follows:
(Folder) Sample name
\n-> (File) (SampleName)_filtered_feature_bc_matrix.h5
-> (Folder) (SampleName)
->-> (File) tissue_lowres_image.png
->-> (File) scalefactors_json.json
->-> (File) tissue_positions_list.csv


special credits to Alfred Kao for providing the code for workflow and many helps
