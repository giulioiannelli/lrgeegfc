# Patient layout normalization plan

## Pat_02
_data/raw/stereoeeg_patients/Pat_02_

| # | kind | src | dst | notes |
|---|------|-----|-----|-------|
| 1 | mkdir | `` | `resting` | create resting/ for canonical phase files |
| 2 | rename | `Resting_PreTask_Data.mat` | `resting/rest_pre.mat` | rest_pre: Resting_PreTask_Data.mat |
| 3 | rename | `Resting_PostTask_Data.mat` | `resting/rest_post.mat` | rest_post: Resting_PostTask_Data.mat |
| 4 | mkdir | `` | `task` | create task/ for canonical phase files |
| 5 | rename | `Task_learning_Data_TimeSeries.mat` | `task/task_learn.mat` | task_learn: Task_learning_Data_TimeSeries.mat |
| 6 | rename | `Task_test_Data_TimeSeries.mat` | `task/task_test.mat` | task_test: Task_test_Data_TimeSeries.mat |
| 7 | mkdir | `` | `implant` | create implant/ |
| 8 | rename | `Implant_pat_2.xlsx` | `implant/implant_pat_02.xlsx` | xlsx: Implant_pat_2.xlsx |
| 9 | generate_csv | `implant/implant_pat_02.xlsx` | `implant_pat_02.csv` | build implant_pat_NN.csv from xlsx |
| 10 | delete | `Implant_pat_02.csv` | `` | legacy implant csv: Implant_pat_02.csv |
| 11 | verify_labels | `channel_labels.txt` | `channel_labels.csv` | verify channel_labels.csv == channel_labels.txt |
| 12 | delete | `channel_labels.txt` | `` | duplicate channel_labels.txt (channel_labels.txt) |
| 13 | delete | `channel_labels.mat` | `` | duplicate channel_labels.mat at patient root |
| 14 | delete | `rsPre.mat` | `` | legacy phase symlink: rsPre.mat |
| 15 | delete | `rsPost.mat` | `` | legacy phase symlink: rsPost.mat |
| 16 | delete | `taskLearn.mat` | `` | legacy phase symlink: taskLearn.mat |
| 17 | delete | `taskTest.mat` | `` | legacy phase symlink: taskTest.mat |
| 18 | write_provenance | `` | `provenance.md` | patient provenance record |


## Pat_03
_data/raw/stereoeeg_patients/Pat_03_

| # | kind | src | dst | notes |
|---|------|-----|-----|-------|
| 1 | mkdir | `` | `resting` | create resting/ for canonical phase files |
| 2 | rename | `RestingPreTask_Data.mat` | `resting/rest_pre.mat` | rest_pre: RestingPreTask_Data.mat |
| 3 | rename | `RestingPostTask_Data.mat` | `resting/rest_post.mat` | rest_post: RestingPostTask_Data.mat |
| 4 | mkdir | `` | `task` | create task/ for canonical phase files |
| 5 | rename | `Task_learning_Data_TimeSeries.mat` | `task/task_learn.mat` | task_learn: Task_learning_Data_TimeSeries.mat |
| 6 | rename | `Task_test_Data_TimeSeries.mat` | `task/task_test.mat` | task_test: Task_test_Data_TimeSeries.mat |
| 7 | mkdir | `` | `implant` | create implant/ |
| 8 | rename | `Implant_pat_3.xlsx` | `implant/implant_pat_03.xlsx` | xlsx: Implant_pat_3.xlsx |
| 9 | generate_csv | `implant/implant_pat_03.xlsx` | `implant_pat_03.csv` | build implant_pat_NN.csv from xlsx |
| 10 | delete | `Implant_pat_03.csv` | `` | legacy implant csv: Implant_pat_03.csv |
| 11 | verify_labels | `channel_labels.txt` | `channel_labels.csv` | verify channel_labels.csv == channel_labels.txt |
| 12 | delete | `channel_labels.txt` | `` | duplicate channel_labels.txt (channel_labels.txt) |
| 13 | delete | `channel_labels.mat` | `` | duplicate channel_labels.mat at patient root |
| 14 | delete | `rsPre.mat` | `` | legacy phase symlink: rsPre.mat |
| 15 | delete | `rsPost.mat` | `` | legacy phase symlink: rsPost.mat |
| 16 | delete | `taskLearn.mat` | `` | legacy phase symlink: taskLearn.mat |
| 17 | delete | `taskTest.mat` | `` | legacy phase symlink: taskTest.mat |
| 18 | write_provenance | `` | `provenance.md` | patient provenance record |


## Pat_05
_data/raw/stereoeeg_patients/Pat_05_

| # | kind | src | dst | notes |
|---|------|-----|-----|-------|
| 1 | mkdir | `` | `resting` | create resting/ for canonical phase files |
| 2 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Pre_Data.mat` | `resting/rest_pre.mat` | rest_pre: PRE_STIM_resting_Pre_Data.mat |
| 3 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Post_Data.mat` | `resting/rest_post.mat` | rest_post: PRE_STIM_resting_Post_Data.mat |
| 4 | mkdir | `` | `task` | create task/ for canonical phase files |
| 5 | rename | `0.53to300Hz_task/PRE_STIM_task_learning_Data.mat` | `task/task_learn.mat` | task_learn: PRE_STIM_task_learning_Data.mat |
| 6 | rename | `0.53to300Hz_task/PRE_STIM_task_test_Data.mat` | `task/task_test.mat` | task_test: PRE_STIM_task_test_Data.mat |
| 7 | mkdir | `` | `implant` | create implant/ |
| 8 | rename | `Implant_locations/Implant_pat_5.xlsx` | `implant/implant_pat_05.xlsx` | xlsx: Implant_pat_5.xlsx |
| 9 | generate_csv | `implant/implant_pat_05.xlsx` | `implant_pat_05.csv` | build implant_pat_NN.csv from xlsx |
| 10 | delete | `Implant_pat_05.csv` | `` | legacy implant csv: Implant_pat_05.csv |
| 11 | verify_labels | `0.53to300Hz_task/channel_labels.txt` | `channel_labels.csv` | verify channel_labels.csv == 0.53to300Hz_task/channel_labels.txt |
| 12 | verify_labels | `0.53to300Hz_resting/channel_labels.txt` | `0.53to300Hz_task/channel_labels.txt` | verify 0.53to300Hz_resting/channel_labels.txt == 0.53to300Hz_task/channel_labels.txt |
| 13 | delete | `0.53to300Hz_task/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_task/channel_labels.txt) |
| 14 | delete | `0.53to300Hz_resting/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_resting/channel_labels.txt) |
| 15 | delete | `rsPre.mat` | `` | legacy phase symlink: rsPre.mat |
| 16 | delete | `rsPost.mat` | `` | legacy phase symlink: rsPost.mat |
| 17 | delete | `taskLearn.mat` | `` | legacy phase symlink: taskLearn.mat |
| 18 | delete | `taskTest.mat` | `` | legacy phase symlink: taskTest.mat |
| 19 | rmdir | `0.53to300Hz_resting` | `` | remove empty vendor subdir: 0.53to300Hz_resting/ |
| 20 | rmdir | `0.53to300Hz_task` | `` | remove empty vendor subdir: 0.53to300Hz_task/ |
| 21 | rmdir | `Implant_locations` | `` | remove empty vendor subdir: Implant_locations/ |
| 22 | write_provenance | `` | `provenance.md` | patient provenance record |


## Pat_06
_data/raw/stereoeeg_patients/Pat_06_

| # | kind | src | dst | notes |
|---|------|-----|-----|-------|
| 1 | mkdir | `` | `resting` | create resting/ for canonical phase files |
| 2 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Pre_Data.mat` | `resting/rest_pre.mat` | rest_pre: PRE_STIM_resting_Pre_Data.mat |
| 3 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Post_Data.mat` | `resting/rest_post.mat` | rest_post: PRE_STIM_resting_Post_Data.mat |
| 4 | mkdir | `` | `task` | create task/ for canonical phase files |
| 5 | rename | `0.53to300Hz_task/PRE_STIM_task_learning_Data.mat` | `task/task_learn.mat` | task_learn: PRE_STIM_task_learning_Data.mat |
| 6 | rename | `0.53to300Hz_task/PRE_STIM_task_test_Data.mat` | `task/task_test.mat` | task_test: PRE_STIM_task_test_Data.mat |
| 7 | mkdir | `` | `implant` | create implant/ |
| 8 | rename | `Implant_locations/Implant_pat_6.xlsx` | `implant/implant_pat_06.xlsx` | xlsx: Implant_pat_6.xlsx |
| 9 | generate_csv | `implant/implant_pat_06.xlsx` | `implant_pat_06.csv` | build implant_pat_NN.csv from xlsx |
| 10 | generate_labels | `0.53to300Hz_task/channel_labels.txt` | `channel_labels.csv` | channel_labels.csv from 0.53to300Hz_task/channel_labels.txt |
| 11 | verify_labels | `0.53to300Hz_resting/channel_labels.txt` | `0.53to300Hz_task/channel_labels.txt` | verify 0.53to300Hz_resting/channel_labels.txt == 0.53to300Hz_task/channel_labels.txt |
| 12 | delete | `0.53to300Hz_task/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_task/channel_labels.txt) |
| 13 | delete | `0.53to300Hz_resting/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_resting/channel_labels.txt) |
| 14 | rmdir | `0.53to300Hz_resting` | `` | remove empty vendor subdir: 0.53to300Hz_resting/ |
| 15 | rmdir | `0.53to300Hz_task` | `` | remove empty vendor subdir: 0.53to300Hz_task/ |
| 16 | rmdir | `Implant_locations` | `` | remove empty vendor subdir: Implant_locations/ |
| 17 | write_provenance | `` | `provenance.md` | patient provenance record |


## Pat_07
_data/raw/stereoeeg_patients/Pat_07_

| # | kind | src | dst | notes |
|---|------|-----|-----|-------|
| 1 | mkdir | `` | `resting` | create resting/ for canonical phase files |
| 2 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Pre_Data.mat` | `resting/rest_pre.mat` | rest_pre: PRE_STIM_resting_Pre_Data.mat |
| 3 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Post_Data.mat` | `resting/rest_post.mat` | rest_post: PRE_STIM_resting_Post_Data.mat |
| 4 | mkdir | `` | `task` | create task/ for canonical phase files |
| 5 | rename | `0.53to300Hz_task/PRE_STIM_task_learning_Data.mat` | `task/task_learn.mat` | task_learn: PRE_STIM_task_learning_Data.mat |
| 6 | rename | `0.53to300Hz_task/PRE_STIM_task_test_Data.mat` | `task/task_test.mat` | task_test: PRE_STIM_task_test_Data.mat |
| 7 | mkdir | `` | `implant` | create implant/ |
| 8 | rename | `Implant_locations/Implant_pat_7.xlsx` | `implant/implant_pat_07.xlsx` | xlsx: Implant_pat_7.xlsx |
| 9 | generate_csv | `implant/implant_pat_07.xlsx` | `implant_pat_07.csv` | build implant_pat_NN.csv from xlsx |
| 10 | delete | `Implant_pat_07.csv` | `` | legacy implant csv: Implant_pat_07.csv |
| 11 | verify_labels | `0.53to300Hz_task/channel_labels.txt` | `channel_labels.csv` | verify channel_labels.csv == 0.53to300Hz_task/channel_labels.txt |
| 12 | verify_labels | `0.53to300Hz_resting/channel_labels.txt` | `0.53to300Hz_task/channel_labels.txt` | verify 0.53to300Hz_resting/channel_labels.txt == 0.53to300Hz_task/channel_labels.txt |
| 13 | delete | `0.53to300Hz_task/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_task/channel_labels.txt) |
| 14 | delete | `0.53to300Hz_resting/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_resting/channel_labels.txt) |
| 15 | delete | `rsPre.mat` | `` | legacy phase symlink: rsPre.mat |
| 16 | delete | `rsPost.mat` | `` | legacy phase symlink: rsPost.mat |
| 17 | delete | `taskLearn.mat` | `` | legacy phase symlink: taskLearn.mat |
| 18 | delete | `taskTest.mat` | `` | legacy phase symlink: taskTest.mat |
| 19 | rmdir | `0.53to300Hz_resting` | `` | remove empty vendor subdir: 0.53to300Hz_resting/ |
| 20 | rmdir | `0.53to300Hz_task` | `` | remove empty vendor subdir: 0.53to300Hz_task/ |
| 21 | rmdir | `Implant_locations` | `` | remove empty vendor subdir: Implant_locations/ |
| 22 | write_provenance | `` | `provenance.md` | patient provenance record |


## Pat_08
_data/raw/stereoeeg_patients/Pat_08_

| # | kind | src | dst | notes |
|---|------|-----|-----|-------|
| 1 | mkdir | `` | `resting` | create resting/ for canonical phase files |
| 2 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Pre_Data.mat` | `resting/rest_pre.mat` | rest_pre: PRE_STIM_resting_Pre_Data.mat |
| 3 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Post_Data.mat` | `resting/rest_post.mat` | rest_post: PRE_STIM_resting_Post_Data.mat |
| 4 | mkdir | `` | `task` | create task/ for canonical phase files |
| 5 | rename | `0.53to300Hz_task/PRE_STIM_task_learning_Data.mat` | `task/task_learn.mat` | task_learn: PRE_STIM_task_learning_Data.mat |
| 6 | rename | `0.53to300Hz_task/PRE_STIM_task_test_Data.mat` | `task/task_test.mat` | task_test: PRE_STIM_task_test_Data.mat |
| 7 | mkdir | `` | `implant` | create implant/ |
| 8 | rename | `Implant_locations/Implant_pat_8.xlsx` | `implant/implant_pat_08.xlsx` | xlsx: Implant_pat_8.xlsx |
| 9 | generate_csv | `implant/implant_pat_08.xlsx` | `implant_pat_08.csv` | build implant_pat_NN.csv from xlsx |
| 10 | delete | `Implant_pat_08.csv` | `` | legacy implant csv: Implant_pat_08.csv |
| 11 | verify_labels | `0.53to300Hz_task/channel_labels.txt` | `channel_labels.csv` | verify channel_labels.csv == 0.53to300Hz_task/channel_labels.txt |
| 12 | verify_labels | `0.53to300Hz_resting/channel_labels.txt` | `0.53to300Hz_task/channel_labels.txt` | verify 0.53to300Hz_resting/channel_labels.txt == 0.53to300Hz_task/channel_labels.txt |
| 13 | delete | `0.53to300Hz_task/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_task/channel_labels.txt) |
| 14 | delete | `0.53to300Hz_resting/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_resting/channel_labels.txt) |
| 15 | delete | `rsPre.mat` | `` | legacy phase symlink: rsPre.mat |
| 16 | delete | `rsPost.mat` | `` | legacy phase symlink: rsPost.mat |
| 17 | delete | `taskLearn.mat` | `` | legacy phase symlink: taskLearn.mat |
| 18 | delete | `taskTest.mat` | `` | legacy phase symlink: taskTest.mat |
| 19 | rmdir | `0.53to300Hz_resting` | `` | remove empty vendor subdir: 0.53to300Hz_resting/ |
| 20 | rmdir | `0.53to300Hz_task` | `` | remove empty vendor subdir: 0.53to300Hz_task/ |
| 21 | rmdir | `Implant_locations` | `` | remove empty vendor subdir: Implant_locations/ |
| 22 | write_provenance | `` | `provenance.md` | patient provenance record |


## Pat_10
_data/raw/stereoeeg_patients/Pat_10_

| # | kind | src | dst | notes |
|---|------|-----|-----|-------|
| 1 | mkdir | `` | `resting` | create resting/ for canonical phase files |
| 2 | rename | `0.53to300Hz_resting/PRE_STIM_restingPre_Data.mat` | `resting/rest_pre.mat` | rest_pre: PRE_STIM_restingPre_Data.mat |
| 3 | rename | `0.53to300Hz_resting/PRE_STIM_restingPost_Data.mat` | `resting/rest_post.mat` | rest_post: PRE_STIM_restingPost_Data.mat |
| 4 | mkdir | `` | `task` | create task/ for canonical phase files |
| 5 | rename | `0.53to300Hz_task/PRE_STIM_task_learning_Data.mat` | `task/task_learn.mat` | task_learn: PRE_STIM_task_learning_Data.mat |
| 6 | rename | `0.53to300Hz_task/PRE_STIM_task_test_Data.mat` | `task/task_test.mat` | task_test: PRE_STIM_task_test_Data.mat |
| 7 | mkdir | `` | `implant` | create implant/ |
| 8 | rename | `Implant_pat_10.xlsx` | `implant/implant_pat_10.xlsx` | xlsx: Implant_pat_10.xlsx |
| 9 | generate_csv | `implant/implant_pat_10.xlsx` | `implant_pat_10.csv` | build implant_pat_NN.csv from xlsx |
| 10 | generate_labels | `0.53to300Hz_task/channel_labels.txt` | `channel_labels.csv` | channel_labels.csv from 0.53to300Hz_task/channel_labels.txt |
| 11 | delete | `0.53to300Hz_task/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_task/channel_labels.txt) |
| 12 | rmdir | `0.53to300Hz_resting` | `` | remove empty vendor subdir: 0.53to300Hz_resting/ |
| 13 | rmdir | `0.53to300Hz_task` | `` | remove empty vendor subdir: 0.53to300Hz_task/ |
| 14 | write_provenance | `` | `provenance.md` | patient provenance record |


## Pat_13
_data/raw/stereoeeg_patients/Pat_13_

| # | kind | src | dst | notes |
|---|------|-----|-----|-------|
| 1 | mkdir | `` | `resting` | create resting/ for canonical phase files |
| 2 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Pre_Data.mat` | `resting/rest_pre.mat` | rest_pre: PRE_STIM_resting_Pre_Data.mat |
| 3 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Post_Data.mat` | `resting/rest_post.mat` | rest_post: PRE_STIM_resting_Post_Data.mat |
| 4 | mkdir | `` | `task` | create task/ for canonical phase files |
| 5 | rename | `0.53to300Hz_task/PRE_STIM_task_learning_Data.mat` | `task/task_learn.mat` | task_learn: PRE_STIM_task_learning_Data.mat |
| 6 | rename | `0.53to300Hz_task/PRE_STIM_task_test_Data.mat` | `task/task_test.mat` | task_test: PRE_STIM_task_test_Data.mat |
| 7 | mkdir | `` | `implant` | create implant/ |
| 8 | rename | `Implant_locations/Implant_pat_13.xlsx` | `implant/implant_pat_13.xlsx` | xlsx: Implant_pat_13.xlsx |
| 9 | generate_csv | `implant/implant_pat_13.xlsx` | `implant_pat_13.csv` | build implant_pat_NN.csv from xlsx |
| 10 | generate_labels | `0.53to300Hz_task/channel_labels.txt` | `channel_labels.csv` | channel_labels.csv from 0.53to300Hz_task/channel_labels.txt |
| 11 | delete | `0.53to300Hz_task/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_task/channel_labels.txt) |
| 12 | rmdir | `0.53to300Hz_resting` | `` | remove empty vendor subdir: 0.53to300Hz_resting/ |
| 13 | rmdir | `0.53to300Hz_task` | `` | remove empty vendor subdir: 0.53to300Hz_task/ |
| 14 | rmdir | `Implant_locations` | `` | remove empty vendor subdir: Implant_locations/ |
| 15 | write_provenance | `` | `provenance.md` | patient provenance record |


## Pat_14
_data/raw/stereoeeg_patients/Pat_14_

| # | kind | src | dst | notes |
|---|------|-----|-----|-------|
| 1 | mkdir | `` | `resting` | create resting/ for canonical phase files |
| 2 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Pre_Data.mat` | `resting/rest_pre.mat` | rest_pre: PRE_STIM_resting_Pre_Data.mat |
| 3 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Post_Data.mat` | `resting/rest_post.mat` | rest_post: PRE_STIM_resting_Post_Data.mat |
| 4 | mkdir | `` | `task` | create task/ for canonical phase files |
| 5 | rename | `0.53to300Hz_task/PRE_STIM_task_learning_Data.mat` | `task/task_learn.mat` | task_learn: PRE_STIM_task_learning_Data.mat |
| 6 | rename | `0.53to300Hz_task/PRE_STIM_task_test_Data.mat` | `task/task_test.mat` | task_test: PRE_STIM_task_test_Data.mat |
| 7 | mkdir | `` | `implant` | create implant/ |
| 8 | rename | `Implant_locations/implant_pat_14.xlsx` | `implant/implant_pat_14.xlsx` | xlsx: implant_pat_14.xlsx |
| 9 | generate_csv | `implant/implant_pat_14.xlsx` | `implant_pat_14.csv` | build implant_pat_NN.csv from xlsx |
| 10 | generate_labels | `0.53to300Hz_task/channel_labels.txt` | `channel_labels.csv` | channel_labels.csv from 0.53to300Hz_task/channel_labels.txt |
| 11 | verify_labels | `0.53to300Hz_resting/channel_labels.txt` | `0.53to300Hz_task/channel_labels.txt` | verify 0.53to300Hz_resting/channel_labels.txt == 0.53to300Hz_task/channel_labels.txt |
| 12 | delete | `0.53to300Hz_task/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_task/channel_labels.txt) |
| 13 | delete | `0.53to300Hz_resting/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_resting/channel_labels.txt) |
| 14 | rmdir | `0.53to300Hz_resting` | `` | remove empty vendor subdir: 0.53to300Hz_resting/ |
| 15 | rmdir | `0.53to300Hz_task` | `` | remove empty vendor subdir: 0.53to300Hz_task/ |
| 16 | rmdir | `Implant_locations` | `` | remove empty vendor subdir: Implant_locations/ |
| 17 | write_provenance | `` | `provenance.md` | patient provenance record |


## Pat_15
_data/raw/stereoeeg_patients/Pat_15_

**Warnings:**
- ⚠ non-canonical implant xlsx stem: found implant_CM.xlsx (expected implant_pat_15.xlsx); renaming it during migration

| # | kind | src | dst | notes |
|---|------|-----|-----|-------|
| 1 | mkdir | `` | `resting` | create resting/ for canonical phase files |
| 2 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Pre_Data.mat` | `resting/rest_pre.mat` | rest_pre: PRE_STIM_resting_Pre_Data.mat |
| 3 | rename | `0.53to300Hz_resting/PRE_STIM_resting_Post_Data.mat` | `resting/rest_post.mat` | rest_post: PRE_STIM_resting_Post_Data.mat |
| 4 | mkdir | `` | `task` | create task/ for canonical phase files |
| 5 | rename | `0.53to300Hz_task/PRE_STIM_task_learning_Data.mat` | `task/task_learn.mat` | task_learn: PRE_STIM_task_learning_Data.mat |
| 6 | rename | `0.53to300Hz_task/PRE_STIM_task_test_Data.mat` | `task/task_test.mat` | task_test: PRE_STIM_task_test_Data.mat |
| 7 | mkdir | `` | `implant` | create implant/ |
| 8 | rename | `Implant_locations/implant_CM.xlsx` | `implant/implant_pat_15.xlsx` | xlsx: implant_CM.xlsx |
| 9 | generate_csv | `implant/implant_pat_15.xlsx` | `implant_pat_15.csv` | build implant_pat_NN.csv from xlsx |
| 10 | generate_labels | `0.53to300Hz_task/channel_labels.txt` | `channel_labels.csv` | channel_labels.csv from 0.53to300Hz_task/channel_labels.txt |
| 11 | verify_labels | `0.53to300Hz_resting/channel_labels.txt` | `0.53to300Hz_task/channel_labels.txt` | verify 0.53to300Hz_resting/channel_labels.txt == 0.53to300Hz_task/channel_labels.txt |
| 12 | delete | `0.53to300Hz_task/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_task/channel_labels.txt) |
| 13 | delete | `0.53to300Hz_resting/channel_labels.txt` | `` | duplicate channel_labels.txt (0.53to300Hz_resting/channel_labels.txt) |
| 14 | rmdir | `0.53to300Hz_resting` | `` | remove empty vendor subdir: 0.53to300Hz_resting/ |
| 15 | rmdir | `0.53to300Hz_task` | `` | remove empty vendor subdir: 0.53to300Hz_task/ |
| 16 | rmdir | `Implant_locations` | `` | remove empty vendor subdir: Implant_locations/ |
| 17 | write_provenance | `` | `provenance.md` | patient provenance record |

