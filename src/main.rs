#[macro_use] extern crate rocket;

mod paste_id;
use std::env;
use paste_id::PasteId;
use sha256::try_digest;
use InSilicoPCR::in_silico_pcr;
use rocket::tokio::fs::File;
use std::path::Path;
use std::path::PathBuf;
use crate::rocket::futures::TryFutureExt;
use rocket::fs::NamedFile;
use anyhow::anyhow;
use rocket::serde::Deserialize;
use rocket::data::Limits;
use rocket::http::uri::Absolute;
use rand::RngCore;
use chacha20poly1305::ChaCha20Poly1305; 
use chacha20poly1305::aead::NewAead;
use chacha20poly1305::aead::stream::EncryptorBE32;
use rocket::{catchers, get, http::Status, http::Method};
use rand::rngs::OsRng;
use rocket::fs::relative;
use rocket::http::Header;
use rocket::response::status::NotFound;
use rocket::tokio::fs::{self};
use rocket::tokio::io::{AsyncReadExt, AsyncWriteExt};
use rocket::routes;
use rocket::fs::{FileServer, TempFile};
use rocket::form::{Form, FromForm};
use rocket::request::FromParam;
use rocket::{Request, Response, response::Redirect};
use reqwest;
use rocket::serde::json::Value;
use rocket::fairing::{Fairing, Info, Kind};
use rocket_governor::{Quota, RocketGovernable, RocketGovernor};
use rocket_governor::rocket_governor_catcher;

pub struct CORS;

#[derive(Debug, Deserialize)]
struct Token {
    challenge_ts: String,
    hostname: String,
    success: bool,
    #[serde(default)]
    error_codes: Vec<String>,
}


#[rocket::get("/")]
fn keyfetch() -> String {
   match env::var("RECAPTCHA_KEY") {
       Ok(secret_key) => format!("{}", secret_key),
       Err(_) => "Secret key not found".to_string(),
       }
}


#[rocket::post("/verify_recaptcha", data="<token>")]
async fn verify_recaptcha(token: String) -> Result<String, Status> {
   let secret_key = keyfetch();
   println!("just for now secret key is {:?}", secret_key);
   if token == "" { return Err(Status::Unauthorized); };
   let url = format!("https://www.google.com/recaptcha/api/siteverify?secret={}&response={}", secret_key, token);
   println!("the url is {:?}", &url);
   let response = reqwest::get(&url).await.map_err(|_| Status::InternalServerError)?;
   match response.json::<Token>().await {
      Ok(resp) => {
             //let body: Value = resp.json().await.map_err(|_| Status::InternalServerError)?;
	     println!("so resp success is {:?}", &resp.success);
	     if resp.success {
	         return Ok("Success".into());
	         }
             else {
	         return Err(Status::Unauthorized);
		 }
             }
      Err(err) => {
             println!("Request failed: {}", err.to_string());
	     return Err(Status::InternalServerError);
	     }
      }
}


#[rocket::async_trait]
impl Fairing for CORS {
  fn info(&self) -> Info {
     Info {
        name: "Add CORS headers to responses",
	kind: Kind::Response,
	}
  }
  async fn on_response<'r>(&self, request: &'r Request<'_>, response: &mut Response<'r>) {
    //let allowed_origins = "http://127.0.0.1:8000";
    let allowed_origins = "http://127.0.0.1:8000";
    if request.method() == Method::Options {
      response.set_status(Status::NoContent);
      response.set_header(Header::new(
         "Access-Control-Allow-Methods",
	 "POST, GET, DELETE, OPTIONS",
	 ));
	 response.set_header(Header::new("Access-Control-Allow-Headers", "Content-Type"));
	 }
	 response.set_header(Header::new(
	    "Access-Control-Allow-Origin",
	    allowed_origins,
	    ));
	 response.set_header(Header::new("Access-Control-Allow-Credentials", "true"));
	 }
}

// In a real application, these would be retrieved dynamically from a config.
const ID_LENGTH: usize = 4;
const HOST: Absolute<'static> = uri!("http://127.0.0.1:8080");


#[derive(FromForm)]
#[derive(Debug)]
struct Pcrrequest<'f> {
    seqfile: TempFile<'f>,
    #[field(default = false)]
    option_m: bool,
    #[field(default = false)]
    option_c: bool,
    #[field(default = 3001)]
    #[field(validate = range(40..3001))]
    option_l: usize,
    primers: TempFile<'f>,
    recaptcha: String,
}

pub struct RateLimitGuard;

impl<'r> RocketGovernable<'r> for RateLimitGuard {
    fn quota(_method: rocket_governor::Method, _route_name: &str) -> Quota {
        Quota::per_second(Self::nonzero(1u32))
    }
}


#[rocket::get("/")]
pub async fn serve() -> Option<NamedFile> {
    let path = Path::new(relative!("assets/index.html"));
    NamedFile::open(path).await.ok()
}

#[get("/assets/<file..>")]
pub async fn servefiles(file: PathBuf) -> Option<NamedFile> {
    println!("so file is {:?}", &file);
    NamedFile::open(Path::new("assets/").join(file)).await.ok()
}

#[get("/route_example")]
fn route_example(_limitguard: RocketGovernor<RateLimitGuard>) -> Status {
    Status::Ok
}

#[rocket::get("/retrieve/<result>/<shasum>")]
async fn retrieve(result: PathBuf, shasum: &str) -> Result<NamedFile, NotFound<String>> {
    let path = Path::new(relative!("/")).join(&result);
    let digest = try_digest(&path).unwrap();
    if digest == shasum {
           NamedFile::open(&path).await.map_err(|e| NotFound(e.to_string()))
       } else {
           NamedFile::open("").await.map_err(|e| NotFound(e.to_string()))
       }
}

pub async fn encrypt_large_file(source_file: PathBuf, dist_file: PathBuf, k3: [u8; 32], n3: [u8; 7]) -> Result<(), anyhow::Error> {
    let aead = ChaCha20Poly1305::new(k3.as_ref().into());
    let mut stream_encryptor = EncryptorBE32::from_aead(aead, n3.as_ref().into());
    const BUFFER_LEN: usize = 500;
    let mut buffer = [0u8; BUFFER_LEN];
    let mut source_file = File::open(source_file).await.expect("source file not open");
    let mut dist_file = File::create(dist_file).await.expect("destination file not created");
    loop {
       let read_count = source_file.read(&mut buffer).await?;
       if read_count == BUFFER_LEN {
          let ciphertext = stream_encryptor
	     .encrypt_next(buffer.as_slice())
	     .map_err(|err| anyhow!("Encrypting large file error {}", err))?;
	  dist_file.write_all(&ciphertext).await?;
       } else {
          let ciphertext = stream_encryptor
	     .encrypt_last(&buffer[..read_count])
             .map_err(|err| anyhow!("Encrypting large file last slice error {}", err))?;
	  dist_file.write_all(&ciphertext).await?;
	  break;
      }
    }
    Ok(())
}

#[rocket::post("/upload", data = "<Pcrrequest>")]
async fn upload(mut Pcrrequest: Form<Pcrrequest<'_>>) -> Result<NamedFile, std::io::Error> {
    let id1 = PasteId::new(ID_LENGTH);
    let id2 = PasteId::new(ID_LENGTH);
    let id3 = PasteId::new(ID_LENGTH);
    let id4 = PasteId::new(ID_LENGTH);
    let id5 = PasteId::new(ID_LENGTH);
    let mut k2 = [0u8; 32];
    OsRng.fill_bytes(&mut k2);
    let mut n2 = [0u8; 7];
    OsRng.fill_bytes(&mut n2);
    let mut k3 = [0u8; 32];
    OsRng.fill_bytes(&mut k3);
    let mut n3 = [0u8; 7];
    OsRng.fill_bytes(&mut n3);
    println!("this is the recaptcha string {:?}", &Pcrrequest.recaptcha);
    let token = &Pcrrequest.recaptcha;
    let verify_token = verify_recaptcha(token.to_string());
    println!("this is the verify token {:?}", &verify_token.await.unwrap());
    Pcrrequest.seqfile.persist_to(id1.file_path()).await?;
    encrypt_large_file(id1.file_path(), id5.file_path(), k2, n2).await.expect("issue with function");
    fs::remove_file(id1.file_path()).await.expect(format!("unable to delete {:?}", id1.file_path().display()).as_str());
    Pcrrequest.primers.persist_to(id2.file_path()).await?;
    encrypt_large_file(id2.file_path(), id4.file_path(), k3, n3).await.expect("issue with second function");
    fs::remove_file(id2.file_path()).await.expect(format!("unable to delete {:?}", id2.file_path().display()).as_str());
    let result = in_silico_pcr(&id5.file_path(), Pcrrequest.option_m, Pcrrequest.option_c, Pcrrequest.option_l, &id4.file_path(), &k3, &n3, &k2, &n2, &id3.file_path());
    let shasum = try_digest(id3.file_path()).unwrap();
    let _finsy = retrieve(id3.file_path(),shasum.as_str());
    fs::remove_file(id4.file_path()).await.expect(format!("unable to delete {:?}", id4.file_path().display()).as_str());
    fs::remove_file(id5.file_path()).await.expect(format!("unable to delete {:?}", id5.file_path().display()).as_str());
    match NamedFile::open(id3.file_path()).await {
      Ok(file) => Ok(file),
      Err(e) => Err(e),
      } 
}


#[launch]
async fn rocket() -> _ {
    let rocket = rocket::build()
       .attach(CORS)
       .mount("/", routes![serve, upload, retrieve, route_example, verify_recaptcha, servefiles])
       .register("/", catchers!(rocket_governor_catcher));
    rocket
}
