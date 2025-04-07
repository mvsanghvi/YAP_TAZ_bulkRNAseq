SELECT top.*
FROM main."top_100" AS top
Left join main."TFlink_Homo_sapiens_browse" as TF on top."Gene Symbol"=TF."Gene Symbol"
WHERE TF.Type LIKE 'transcription factor%'
