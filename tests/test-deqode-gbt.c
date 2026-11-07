/* Offline regression tests for deqode GBT integer decoding.
 * No mining threads, sockets, RPC client or submission code is invoked.
 */
#include <limits.h>
#include <stdarg.h>

/* Include the real translation unit so its static decoder is tested directly. */
#include "../algo/deqode/deqode-gate.c"

bool opt_debug = false;

/* Logging is outside this test's scope. */
void applog( int priority, const char *format, ... )
{
   (void) priority;
   (void) format;
}

/* The fixture's acceptance target is 32 zero bytes. Target parsing itself is
 * outside these integer-decoding tests; this stub is linked only into the test.
 */
bool hex2bin( unsigned char *out, const char *hex, const size_t len )
{
   if ( strlen( hex ) != 2 * len ) return false;
   for ( size_t i = 0; i < 2 * len; ++i )
      if ( hex[i] != '0' ) return false;
   memset( out, 0, len );
   return true;
}

static unsigned checks = 0, failures = 0;

static void check( const char *name, bool passed )
{
   ++checks;
   if ( !passed ) ++failures;
   printf( "%s %u - %s\n", passed ? "ok" : "not ok", checks, name );
}

static void integer_case( const char *name, json_t *value,
                          bool accepted, int expected )
{
   json_t *object = json_object();
   if ( value ) json_object_set_new( object, "value", value );
   int out = 12345;
   const bool result = json_int_field( object, "value", &out );
   check( name, result == accepted && out == ( accepted ? expected : 12345 ) );
   json_decref( object );
}

static json_t *valid_template( void )
{
   json_t *sizing = json_pack( "{s:i,s:i,s:i,s:i,s:i,s:s}",
      "n", 128, "r", 64, "w", 4, "gate_multiplier", 24,
      "expected_witnesses", 1,
      "acceptance_target",
      "0000000000000000000000000000000000000000000000000000000000000000" );
   json_t *task = json_pack( "{s:i,s:i,s:o}",
      "id", DEQODE_TASK_ID, "data_bytes", 32, "deqode", sizing );
   return json_pack( "{s:o}", "task", task );
}

static bool decode( const json_t *value, struct work *out )
{
   memset( out, 0, sizeof( *out ) );
   out->data[0] = bswap_32( UINT32_C(0x40100000) );
   return deqode_gbt_decode_task( value, out );
}

static void field_case( const char *name, bool sizing_field,
                        const char *field, json_t *replacement, bool accepted )
{
   json_t *value = valid_template();
   json_t *object = json_object_get( value, "task" );
   if ( sizing_field ) object = json_object_get( object, "deqode" );
   if ( replacement ) json_object_set_new( object, field, replacement );
   else json_object_del( object, field );
   struct work work;
   check( name, decode( value, &work ) == accepted );
   json_decref( value );
}

int main( void )
{
   /* All currently supported miner platforms use a 32-bit int. */
   _Static_assert( INT_MAX == 2147483647, "These wraparound fixtures require 32-bit int" );
   _Static_assert( sizeof( json_int_t ) > sizeof( int ), "Wide Jansson integers required" );

   integer_case( "integer zero accepted", json_integer( 0 ), true, 0 );
   integer_case( "negative in-range integer accepted", json_integer( -17 ), true, -17 );
   integer_case( "INT_MIN accepted", json_integer( INT_MIN ), true, INT_MIN );
   integer_case( "INT_MAX accepted", json_integer( INT_MAX ), true, INT_MAX );
   integer_case( "missing integer rejected without modifying output", NULL, false, 0 );
   integer_case( "null integer rejected without modifying output", json_null(), false, 0 );
   integer_case( "string integer rejected without modifying output", json_string( "1" ), false, 0 );
   integer_case( "real integer rejected without modifying output", json_real( 1.0 ), false, 0 );
   integer_case( "boolean integer rejected without modifying output", json_true(), false, 0 );
   integer_case( "INT_MAX plus one rejected without modifying output",
      json_integer( (json_int_t) INT_MAX + 1 ), false, 0 );
   integer_case( "INT_MIN minus one rejected without modifying output",
      json_integer( (json_int_t) INT_MIN - 1 ), false, 0 );
   integer_case( "large wrapped integer rejected without modifying output",
      json_integer( (json_int_t) 4294967297LL ), false, 0 );

   json_t *value = valid_template();
   struct work work;
   const bool decoded = decode( value, &work );
   check( "valid task and sizing preserved", decoded && work.task_params[0] == 128
      && work.task_params[1] == 64 && work.task_params[2] == 4 && work.task_params[3] == 24 );
   json_decref( value );

   field_case( "wrapped task id rejected", false, "id", json_integer( 4294967297LL ), false );
   const char *fields[] = { "n", "r", "w", "gate_multiplier", "expected_witnesses" };
   const int originals[] = { 128, 64, 4, 24, 1 };
   for ( unsigned i = 0; i < sizeof( fields ) / sizeof( fields[0] ); ++i )
   {
      char name[100];
      snprintf( name, sizeof(name), "wrapped sizing %s rejected", fields[i] );
      field_case( name, true, fields[i],
         json_integer( (json_int_t) 4294967296LL + originals[i] ), false );
   }
   field_case( "missing optional data_bytes accepted", false, "data_bytes", NULL, true );
   field_case( "correct optional data_bytes accepted", false, "data_bytes", json_integer( 32 ), true );
   field_case( "wrong optional data_bytes rejected", false, "data_bytes", json_integer( 31 ), false );
   field_case( "negative optional data_bytes rejected", false, "data_bytes", json_integer( -1 ), false );
   field_case( "wrapped optional data_bytes rejected", false, "data_bytes", json_integer( 4294967328LL ), false );
   field_case( "out-of-range optional data_bytes rejected", false, "data_bytes",
      json_integer( (json_int_t) INT_MAX + 1 ), false );
   field_case( "null optional data_bytes rejected", false, "data_bytes", json_null(), false );
   field_case( "string optional data_bytes rejected", false, "data_bytes", json_string( "32" ), false );
   field_case( "real optional data_bytes rejected", false, "data_bytes", json_real( 32.0 ), false );
   field_case( "boolean optional data_bytes rejected", false, "data_bytes", json_true(), false );
   field_case( "missing required task id rejected", false, "id", NULL, false );
   field_case( "missing required sizing rejected", true, "n", NULL, false );
   field_case( "null required sizing rejected", true, "n", json_null(), false );
   field_case( "string required sizing rejected", true, "n", json_string( "128" ), false );
   field_case( "real required sizing rejected", true, "n", json_real( 128.0 ), false );
   field_case( "invalid dimensions remain rejected", true, "r", json_integer( 129 ), false );
   printf( "1..%u\n", checks );
   printf( "# %u checks, %u failures\n", checks, failures );
   return failures ? 1 : 0;
}
