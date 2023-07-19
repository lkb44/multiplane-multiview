import pandas as pd
import random

# Set the seed for reproducibility
random.seed(42)

# File path of the spreadsheet
file_path = '/Users/leobao/Documents/MultiPlanePipeline/multiplane-multiview/T0-2vT3-4cohort.csv'

# Read the spreadsheet into a pandas DataFrame
df = pd.read_csv(file_path)

# Calculate the number of patients for each cohort
total_patients = len(df)
train_size = int(0.8 * total_patients)
test_size = total_patients - train_size

# Separate responders and non-responders
responders = df[df['Label'] == 1]
non_responders = df[df['Label'] == 0]

# Shuffle the responders and non-responders
responders = responders.sample(frac=1, random_state=42)
non_responders = non_responders.sample(frac=1, random_state=42)

# Assign cohorts to responders
train_responders = responders[:int(train_size * len(responders) / total_patients)]
test_responders = responders[int(train_size * len(responders) / total_patients):]

# Assign cohorts to non-responders
train_non_responders = non_responders[:int(train_size * len(non_responders) / total_patients)]
test_non_responders = non_responders[int(train_size * len(non_responders) / total_patients):]

# Concatenate the training and testing cohorts
train_cohort = pd.concat([train_responders, train_non_responders])
test_cohort = pd.concat([test_responders, test_non_responders])

# Shuffle the combined cohorts
train_cohort = train_cohort.sample(frac=1, random_state=42)
test_cohort = test_cohort.sample(frac=1, random_state=42)

# Update the assigned cohort column in the DataFrame
df.loc[train_cohort.index, 'Cohort'] = 'Train'
df.loc[test_cohort.index, 'Cohort'] = 'Test'

# Save the updated DataFrame to a new CSV file
output_file_path = '/Users/leobao/Documents/MultiPlanePipeline/multiplane-multiview/T0-2vT3-4cohort_split.csv'
df.to_csv(output_file_path, index=False)

print('Cohort split completed. The updated DataFrame has been saved to', output_file_path)