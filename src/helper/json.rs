pub trait FromJsonString: Sized {
    fn from_json_string(json_str: &str) -> Result<Self, serde_json::Error>;
}

impl<T> FromJsonString for T
where
    T: serde::de::DeserializeOwned,
{
    fn from_json_string(json_str: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json_str)
    }
}
